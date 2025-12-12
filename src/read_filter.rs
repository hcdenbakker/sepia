use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use niffler;
use std;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::process;

pub fn tab_to_map(
    filename: String,
    query: &str,
    ratio: f64,
    exact: bool,
) -> (
    std::collections::HashMap<std::string::String, String>,
    std::collections::HashMap<std::string::String, String>,
) {
    let mut class_map = HashMap::new();
    let mut ratio_map = HashMap::new();
    //let f = File::open(filename).expect("classification file not found");
    let (f, _compression) = niffler::from_path(filename).expect("classification file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        let h: Vec<&str> = v[0].split(' ').collect();
        if (v[2].parse::<f64>().unwrap() / v[3].parse::<f64>().unwrap()) >= ratio {
            ratio_map.insert(String::from(h[0]), String::from(v[1]));
        }
        if exact == true {
            if v[1] == query {
                class_map.insert(String::from(h[0]), String::from(v[1]));
            }
        } else {
            if v[1].contains(query)
            /* && (v[4] == accept)*/
            {
                class_map.insert(String::from(h[0]), String::from(v[1]));
            }
        }
    }
    (class_map, ratio_map)
}

pub fn tab_to_set(
    filename: String,
    query: &str,
    ratio: f64,
    exact: bool,
) -> (
    std::collections::HashSet<std::string::String>,
    std::collections::HashSet<std::string::String>,
) {
    let mut class_set = HashSet::new();
    let mut ratio_set = HashSet::new();
    //let f = File::open(filename).expect("classification file not found");
    let (f, _compression) = niffler::from_path(filename).expect("classification file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        let h: Vec<&str> = v[0].split(' ').collect();
        if (v[2].parse::<f64>().unwrap() / v[3].parse::<f64>().unwrap()) >= ratio {
            ratio_set.insert(String::from(h[0]));
        }
        if exact == true {
            if v[1] == query {
                class_set.insert(String::from(h[0]));
            }
        } else {
            if v[1].contains(query)
            /* && (v[4] == accept)*/
            {
                class_set.insert(String::from(h[0]));
            }
        }
    }
    (class_set, ratio_set)
}

pub fn fastq_or_fasta(file_name: &str) -> String {
    let (d, _format) = niffler::from_path(file_name).expect("nibbler choked");
    let mut first_line = String::new();
    io::BufReader::new(d)
        .read_line(&mut first_line)
        .expect("Unable to read line");
    let first_char = first_line.chars().next().unwrap();
    if first_char == '>' {
        "fasta".to_string()
    } else if first_char == '@' {
        "fastq".to_string()
    } else {
        "unknown".to_string()
    }
}

#[allow(unused_assignments)]
pub fn read_filter_pe(
    class_map: std::collections::HashMap<std::string::String, String>,
    ratio_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let f2 = File::open(&filenames[1]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let d2 = MultiGzDecoder::new(f2);
    let iter1 = io::BufReader::new(d1).lines();
    let mut iter2 = io::BufReader::new(d2).lines();
    let mut header1 = "".to_string();
    let mut header2 = "".to_string();
    let mut seq1 = "".to_string();
    let mut seq2 = "".to_string();
    let mut qual1 = "".to_string();
    let mut qual2 = "".to_string();
    let mut excluded = 0;
    let mut included = 0;
    let fq1 = File::create(format!("{}_R1.fq.gz", prefix)).expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    let fq2 = File::create(format!("{}_R2.fq.gz", prefix)).expect("could not create R2!");
    let mut gz2 = GzEncoder::new(fq2, Compression::default());
    for line in iter1 {
        let l = line.unwrap();
        let line2 = iter2.next();
        //maybe check here if headers come from same pair
        if line_count % 4 == 1 {
            match line2 {
                Some(h2) => {
                    header1 = l.to_string();
                    header2 = h2.unwrap().to_string();
                }
                None => break,
            };
        } else if line_count % 4 == 2 {
            match line2 {
                Some(l2) => {
                    seq1 = l.to_owned();
                    seq2 = l2.unwrap().to_owned();
                }
                None => break,
            };
        } else if line_count % 4 == 0 {
            match line2 {
                Some(l2) => {
                    qual1 = l.to_owned();
                    qual2 = l2.unwrap().to_owned();
                    let v: Vec<&str> = header1.split(' ').collect();
                    if exclude == true {
                        if !class_map.contains_key(&v[0][1..]) {
                            gz1.write_all(
                                format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes(),
                            )
                            .expect("could not write R1!");
                            gz2.write_all(
                                format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes(),
                            )
                            .expect("could not write R2!");
                            excluded += 1;
                        }
                    } else {
                        if (class_map.contains_key(&v[0][1..]) && ratio_map.is_empty())
                            || (class_map.is_empty() && ratio_map.contains_key(&v[0][1..]))
                            || (class_map.contains_key(&v[0][1..])
                                && ratio_map.contains_key(&v[0][1..]))
                        {
                            gz1.write_all(
                                format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes(),
                            )
                            .expect("could not write R1!");
                            gz2.write_all(
                                format!("{}\n{}\n+\n{}\n", header2, seq2, qual2).as_bytes(),
                            )
                            .expect("could not write R2!");
                            included += 1;
                        }
                    }
                }
                None => break,
            };
        }
        line_count += 1;
    }
    gz1.finish().expect("Could not close new R1 file");
    gz2.finish().expect("Could not close new R2 file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

//quicker to adjust the function to se
#[allow(unused_assignments)]
pub fn read_filter_se(
    class_map: std::collections::HashMap<std::string::String, String>,
    ratio_map: std::collections::HashMap<std::string::String, String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let mut line_count = 1;
    let f = File::open(&filenames[0]).expect("file not found");
    let d1 = MultiGzDecoder::new(f);
    let iter1 = io::BufReader::new(d1).lines();
    let mut header1 = "".to_string();
    let mut seq1 = "".to_string();
    let mut qual1 = "".to_string();
    let mut excluded = 0;
    let mut included = 0;
    let fq1 = File::create(format!("{}.fq.gz", prefix)).expect("could not create R1!");
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    for line in iter1 {
        let l = line.unwrap();
        if line_count % 4 == 1 {
            header1 = l.to_string();
        } else if line_count % 4 == 2 {
            seq1 = l.to_owned();
        } else if line_count % 4 == 0 {
            qual1 = l.to_owned();
            let v: Vec<&str> = header1.split(' ').collect();
            if exclude == true {
                if !class_map.contains_key(&v[0][1..]) {
                    gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                        .expect("Could not write forward read(-s) to file");
                    excluded += 1;
                }
            } else {
                if (class_map.contains_key(&v[0][1..]) && ratio_map.is_empty())
                    || (class_map.is_empty() && ratio_map.contains_key(&v[0][1..]))
                    || (class_map.contains_key(&v[0][1..]) && ratio_map.contains_key(&v[0][1..]))
                {
                    gz1.write_all(format!("{}\n{}\n+\n{}\n", header1, seq1, qual1).as_bytes())
                        .expect("Could not write reverse read(-s) to file");
                    included += 1;
                }
            }
        }
        line_count += 1;
    }
    gz1.finish().expect("Could not close new read file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            excluded, query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

//needletail version
#[allow(unused_assignments)]
pub fn read_filter_se_nt(
    class_set: std::collections::HashSet<std::string::String>,
    ratio_set: std::collections::HashSet<std::string::String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    let mut intersection_size: usize = 0;
    for _e in class_set.intersection(&ratio_set){
        intersection_size += 1;
    }
    if (class_set.is_empty() && ratio_set.is_empty() && !exclude) || (!class_set.is_empty() && !ratio_set.is_empty() && intersection_size == 0 && !exclude){
        println!("No read-pairs matching the inclusion criteria for {}, exiting with 0 now.", filenames[0]);
        std::process::exit(0x0100);

    }
    let mut total = 0;
    let mut excluded = 0;
    let mut included = 0;
    let format = fastq_or_fasta(filenames[0]);
    println!("{} detected!", format);
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let fq1 = if format == "fastq" {
        File::create(format!("{}.fq.gz", prefix)).expect("could not create output file!")
    } else {
        File::create(format!("{}.fa.gz", prefix)).expect("could not create output file!")
    };
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    while let Some(record1) = reader1.next() {
        total += 1;
        let seqrec1 = record1.expect("invalid record in forward file");
        let v: Vec<&str> = std::str::from_utf8(seqrec1.id())
            .unwrap()
            .split(' ')
            .collect();
        //let v: Vec<&str> = seqrec1.id().split(' ').collect();
        if exclude == true {
            if (!class_set.contains(v[0]) && ratio_set.is_empty())
                || !(class_set.contains(v[0]) && ratio_set.contains(v[0]))
            {
                match seqrec1.format() {
                    needletail::parser::Format::Fasta => gz1
                        .write_all(
                            format!(
                                ">{}\n{}\n",
                                std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string()
                            )
                            .as_bytes(),
                        )
                        .expect("Could not write reverse read(-s) to file"),
                    //} else {
                    needletail::parser::Format::Fastq => gz1
                        .write_all(
                            format!(
                                "@{}\n{}\n+\n{}\n",
                                std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.qual().unwrap())
                                    .unwrap()
                                    .to_string()
                            )
                            .as_bytes(),
                        )
                        .expect("Could not write reverse read(-s) to file"),
                }
                //record1.unwrap().write(&gz1, Some(record1.unwrap().line_ending()));
                excluded += 1;
            }
        } else {
            if (class_set.contains(v[0]) && ratio_set.is_empty())
                || (class_set.is_empty() && ratio_set.contains(v[0]))
                || (class_set.contains(v[0]) && ratio_set.contains(v[0]))
            {
                match seqrec1.format() {
                    needletail::parser::Format::Fasta => gz1
                        .write_all(
                            format!(
                                ">{}\n{}\n",
                                std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string()
                            )
                            .as_bytes(),
                        )
                        .expect("Could not write reverse read(-s) to file"),
                    //} else {
                    needletail::parser::Format::Fastq => gz1
                        .write_all(
                            format!(
                                "@{}\n{}\n+\n{}\n",
                                std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string(),
                                std::str::from_utf8(&seqrec1.qual().unwrap())
                                    .unwrap()
                                    .to_string()
                            )
                            .as_bytes(),
                        )
                        .expect("Could not write reverse read(-s) to file"),
                }
                included += 1;
            }
        }
    }
    gz1.finish().expect("Could not close new read file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            total - excluded,
            query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

pub fn read_filter_pe_nt(
    class_set: std::collections::HashSet<std::string::String>,
    ratio_set: std::collections::HashSet<std::string::String>,
    filenames: Vec<&str>,
    query: &str,
    prefix: &str,
    exclude: bool,
) {
    //if none of the inclusion criteria are met, print there is nothing to include and exit with 0
    //this does not work if a minimum kmer-ratio is needed
    let mut intersection_size: usize = 0;
    for _e in class_set.intersection(&ratio_set){
        intersection_size += 1;
    }
    if (class_set.is_empty() && ratio_set.is_empty() && !exclude) || (!class_set.is_empty() && !ratio_set.is_empty() && intersection_size == 0 && !exclude){
        println!("No read-pairs matching the inclusion criteria for {}, exiting with 0 now.", filenames[0]);
        std::process::exit(0x0100);

    }
    let mut total = 0;
    let mut excluded = 0;
    let mut included = 0;
    //detect format; we determine if we have fasta or fastq, or throw an errer if we cannot
    //recognize the format, or if we have a mix of fasta and fastq
    let format = fastq_or_fasta(filenames[0]);
    let mut formats: HashSet<String> = HashSet::new();
    for f in &filenames {
        formats.insert(fastq_or_fasta(f));
    }
    if formats.len() > 1 || (formats.get(&"unknown".to_string()) == Some(&"unknown".to_string())) {
        eprintln!("sequence format not recognized or mixed fasta and fastq data");
        process::abort();
    }
    println!("{} detected!", format);
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut reader2 = parse_fastx_file(&filenames[1]).expect("invalid path/file");
    let fq1 = if format == "fastq" {
        File::create(format!("{}_R1.fq.gz", prefix)).expect("could not create R1 file!")
    } else {
        File::create(format!("{}_R1.fa.gz", prefix)).expect("could not create R1 file!")
    };
    let mut gz1 = GzEncoder::new(fq1, Compression::default());
    let fq2 = if format == "fastq" {
        File::create(format!("{}_R2.fq.gz", prefix)).expect("could not create R2 file!")
    } else {
        File::create(format!("{}_R2.fa.gz", prefix)).expect("could not create R2 file!")
    };
    let mut gz2 = GzEncoder::new(fq2, Compression::default());
    while let Some(record1) = reader1.next() {
        total += 1;
        let seqrec1 = record1.expect("invalid record in forward file");
        if let Some(record2) = reader2.next() {
            let seqrec2 = record2.expect("invalid record in forward file");
            let v: Vec<&str> = std::str::from_utf8(seqrec1.id())
                .unwrap()
                .split(' ')
                .collect();
            if exclude == true {
                if (!class_set.contains(v[0]) && ratio_set.is_empty())
                    || !(class_set.contains(v[0]) && ratio_set.contains(v[0]))
                {
                    match seqrec1.format() {
                        needletail::parser::Format::Fasta => {
                            gz1.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write forward read(-s) to file");
                            gz2.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    std::str::from_utf8(&seqrec2.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.raw_seq()).unwrap().to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write reverse read(-s) to file")
                        }
                        needletail::parser::Format::Fastq => {
                            gz1.write_all(
                                format!(
                                    "@{}\n{}\n+\n{}\n",
                                    std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.qual().unwrap())
                                        .unwrap()
                                        .to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write forward read(-s) to file");
                            gz2.write_all(
                                format!(
                                    "@{}\n{}\n+\n{}\n",
                                    std::str::from_utf8(&seqrec2.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.raw_seq()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.qual().unwrap())
                                        .unwrap()
                                        .to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write reverse read(-s) to file")
                        }
                    }
                    //record1.unwrap().write(&gz1, Some(record1.unwrap().line_ending()));
                    excluded += 1;
                }
            } else {
                if (class_set.contains(v[0]) && ratio_set.is_empty())
                    || (class_set.is_empty() && ratio_set.contains(v[0]))
                    || (class_set.contains(v[0]) && ratio_set.contains(v[0]))
                {
                    match seqrec1.format() {
                        needletail::parser::Format::Fasta => {
                            gz1.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write forward read(-s) to file");
                            gz2.write_all(
                                format!(
                                    ">{}\n{}\n",
                                    std::str::from_utf8(&seqrec2.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.raw_seq()).unwrap().to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write reverse read(-s) to file")
                        }
                        needletail::parser::Format::Fastq => {
                            gz1.write_all(
                                format!(
                                    "@{}\n{}\n+\n{}\n",
                                    std::str::from_utf8(&seqrec1.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.raw_seq()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec1.qual().unwrap())
                                        .unwrap()
                                        .to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write forward read(-s) to file");
                            gz2.write_all(
                                format!(
                                    "@{}\n{}\n+\n{}\n",
                                    std::str::from_utf8(&seqrec2.id()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.raw_seq()).unwrap().to_string(),
                                    std::str::from_utf8(&seqrec2.qual().unwrap())
                                        .unwrap()
                                        .to_string()
                                )
                                .as_bytes(),
                            )
                            .expect("Could not write reverse read(-s) to file")
                        }
                    }

                    included += 1;
                }
            }
        }
    }
    gz1.finish().expect("Could not close new R1 file");
    gz2.finish().expect("Could not close new R2 file");
    if exclude == true {
        eprintln!(
            "Excluded {} read pairs  with classification containing '{}' from output files",
            total - excluded,
            query
        );
    } else {
        eprintln!(
            "Wrote {} read-pairs with classification containing '{}' to output files",
            included, query
        );
    }
}

//use needletail to make fasta/fastq and compression agnostic
