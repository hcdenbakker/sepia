use super::bit_magic::get_sepia;
use super::kmer;
use needletail::parse_fastx_file;
use niffler;

use flate2::write::GzEncoder;
use flate2::Compression;

use hashbrown::HashMap;
use itertools::Itertools;
use rayon::prelude::*;
use std;

use std::collections::VecDeque;

use std::fs::File;

use std::io;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::str;
use std::sync::Arc;
use std::time::SystemTime;

use bstr::io::BufReadExt;
use bstr::ByteVec;
use sdset::Set;

use hyperloglog::HyperLogLog;
use std::cmp::min;

const TOGGLE: u64 = 0xe37e28c4271b5a2d;

#[inline]
fn sliding_window_minimizers_sepia_vec(
    seq: &str,
    k: usize,
    m: usize,
    b: u32,
    db: &[u32],
) -> (Vec<(u32, usize)>, usize) {
    let capacity_db = db.len();
    let mut report: HashMap<u32, usize> = HashMap::default();
    let mut observations = 0;
    let mut window: VecDeque<(u64, usize)> = VecDeque::new(); //position minimizer
    let mut taxon = 0;
    let mut current_minimizer: u64 = 0;
    let length = seq.len();
    let mut counter = 0;
    let mut i = 1;
    let mut candidate: u64 = 0;
    let mut mask: u64 = 1;
    mask <<= m * 2;
    mask -= 1;
    let toggle = TOGGLE & mask;
    for n in seq.bytes() {
        if kmer::nuc_to_number(n) < 4 {
            candidate <<= 2;
            candidate |= kmer::nuc_to_number(n);
            counter += 1;
            if counter >= m {
                candidate &= mask;
                let mut minimizer = kmer::canonical(candidate, m);
                minimizer ^= toggle;
                while !window.is_empty() && window.back().unwrap().0 > minimizer {
                    window.pop_back(); // we pop the last one
                }
                window.push_back((minimizer, i)); // and make add a pair with the new value at the end
                while window.front().unwrap().1 as isize <= i as isize - k as isize + m as isize - 1
                {
                    window.pop_front(); // pop the first one
                }
                if i >= k {
                    if current_minimizer == window.front().unwrap().0 ^ toggle {
                        //    continue
                        let count = report.entry(taxon).or_insert(0);
                        *count += 1;
                    } else {
                        //continue
                        taxon =
                            get_sepia(&(window.front().unwrap().0 ^ toggle), &db, capacity_db, b);
                        let count = report.entry(taxon).or_insert(0);
                        *count += 1;
                        current_minimizer = window.front().unwrap().0;
                    }
                    observations += 1;
                }
            }
        } else {
            // we have an ambiguous character: reset counter, empty queue
            counter = 0;
            i = 0;
            candidate = 0;
            window.clear();
        }
        if i == length {
            break;
        }
        i += 1;
    }
    let count_vec: Vec<(u32, usize)> = report.into_iter().collect();
    (count_vec, observations)
}

#[inline]
fn sliding_window_minimizers_sepia_prehll(
    seq: &str,
    k: usize,
    m: usize,
    b: u32,
    db: &[u32],
) -> (Vec<(u32, usize)>, usize, Vec<(u32, u64)>) {
    let capacity_db = db.len();
    let l_r = kmer::revcomp(seq);
    //let mut all_kmers = Vec::new();
    let mut report: HashMap<u32, usize> = HashMap::default();
    let mut report_kmers: Vec<(u32, u64)> = Vec::new();
    let mut observations = 0;
    let mut window: VecDeque<(u64, usize)> = VecDeque::new(); //position minimizer
    let mut taxon = 0;
    let mut current_minimizer: u64 = 0;
    let length = seq.len();
    let mut counter = 0;
    let mut i = 1; //counter for minimizer, resets when ambiguous character is enountered
    let mut j = 1; //total seq counter
    let mut candidate: u64 = 0;
    let mut mask: u64 = 1;
    //let mut canonical_kmer = min(&seq[0..31], &l_r[length - (31)..length - 31]);
    mask <<= m * 2;
    mask -= 1;
    let toggle = TOGGLE & mask;
    for n in seq.bytes() {
        //for j in seq[i..i + m].bytes(){
        //let new_char = super::kmer::nuc_to_number(&n);
        /*if j >= k {
            let index = j - k;
            canonical_kmer = min(
                &seq[index..index + k],
                &l_r[length - (index + k)..length - index],
            );
            //println!("{}", canonical_kmer);
        }*/
        if kmer::nuc_to_number(n) < 4 {
            candidate <<= 2;
            candidate |= kmer::nuc_to_number(n);
            counter += 1;
            if counter >= m {
                candidate &= mask;
                //candidate &= seed;
                let mut minimizer = kmer::canonical(candidate, m);
                minimizer ^= toggle;
                while !window.is_empty() && window.back().unwrap().0 > minimizer {
                    window.pop_back(); // we pop the last one
                }
                window.push_back((minimizer, i)); // and make add a pair with the new value at the end
                while window.front().unwrap().1 as isize <= i as isize - k as isize + m as isize - 1
                {
                    window.pop_front(); // pop the first one
                }
                if i >= k {
                    //println!("{}", canonical_kmer);
                    let index = j - k;
                    let canonical_kmer = min(
                        &seq[index..index + k],
                        &l_r[length - (index + k)..length - index],
                    );
                    if current_minimizer == window.front().unwrap().0 {
                        let count = report.entry(taxon).or_insert(0);
                        *count += 1;
                        report_kmers.push((taxon, seahash::hash(canonical_kmer.as_bytes())));
                    } else {
                        taxon =
                            get_sepia(&(window.front().unwrap().0 ^ toggle), &db, capacity_db, b);
                        let count = report.entry(taxon).or_insert(0);
                        *count += 1;
                        report_kmers.push((taxon, seahash::hash(canonical_kmer.as_bytes())));
                        current_minimizer = window.front().unwrap().0;
                    }
                    observations += 1;
                }
            }
        } else {
            // we have an ambiguous character: reset counter, empty queue
            counter = 0;
            i = 0;
            candidate = 0;
            window.clear();
        }
        // }
        //if i == length {
        //    break;
        //}
        i += 1;
        j += 1;
        /*if j == length {
            break;
        }*/
    }
    let count_vec: Vec<(u32, usize)> = report.into_iter().collect();
    (count_vec, observations, report_kmers)
}

#[inline]
pub fn search_index_lca(map: &[String], db: &[u32], bloom_size: u64) -> HashMap<u32, usize> {
    let mut report = HashMap::default();
    for k in map {
        let position = seahash::hash(&k.as_bytes()) % bloom_size;
        let count = report.entry(db[position as usize]).or_insert(0);
        *count += 1;
    }
    report
}

//this does not seem very efficient O(n squared) at it its worst
#[inline]
pub fn score_lineages(
    report: &HashMap<u32, usize>,
    taxonomy: &HashMap<u32, String>,
) -> HashMap<u32, usize> {
    let mut lineage_scores = HashMap::default();
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report.keys() {
        if *key == (0 as u32) || *key == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[&key];
            for key2 in report.keys() {
                if *key2 == (0 as u32) || *key2 == unclassified {
                    continue;
                } else {
                    let lineage2 = &taxonomy[&key2];
                    if lineage1 == lineage2 {
                        let count = lineage_scores.entry(*key).or_insert(0);
                        *count += &report[&key];
                    } else {
                        if lineage1.contains(lineage2) {
                            let count = lineage_scores.entry(*key).or_insert(0);
                            *count += &report[&key2];
                        }
                    }
                }
            }
        }
    }
    lineage_scores
}

#[inline]
pub fn poll_lineages(
    report: &hashbrown::HashMap<u32, usize>,
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
) -> (u32, usize, usize) {
    let mut max_score = 0;
    let mut max_taxon: u32 = 0;
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report.keys() {
        let mut score = 0;
        if *key == (0 as u32) || *key == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[&key];
            for key2 in report.keys() {
                if *key2 == (0 as u32) || *key2 == unclassified {
                    continue;
                } else {
                    let lineage2 = &taxonomy[&key2];
                    //if lineage1 == lineage2 {
                    //    score += &report[&key];
                    //} else {
                    if lineage1.contains(lineage2) {
                        score += &report[&key2];
                    }
                    //}
                }
            }
        }
        if score > max_score {
            max_taxon = *key;
            max_score = score;
        }
        if score == max_score {
            max_taxon =
                super::taxonomy_u32::find_lca_u32(*key, max_taxon, lineage_graph, unclassified);
        }
    }
    (max_taxon, max_score, kmer_length)
}

#[inline]
pub fn poll_lineages_vec(
    report: &Vec<(u32, usize)>,
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
) -> (u32, usize, usize) {
    let mut max_score = 0;
    let mut max_taxon: u32 = 0;
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report {
        let mut score = 0;
        if key.0 == (0 as u32) || key.0 == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[&key.0];
            for key2 in report {
                if key2.0 == (0 as u32) || key2.0 == unclassified {
                    continue;
                } else {
                    let lineage2 = &taxonomy[&key2.0];
                    //if lineage1 == lineage2 {
                    //    score += &report[&key];
                    //} else {
                    if lineage1.contains(lineage2) {
                        score += key2.1;
                    }
                    //}
                }
            }
        }
        if score > max_score {
            max_taxon = key.0;
            max_score = score;
        }
        if score == max_score {
            max_taxon =
                super::taxonomy_u32::find_lca_u32(key.0, max_taxon, lineage_graph, unclassified);
        }
    }
    (max_taxon, max_score, kmer_length)
}

#[inline]
pub fn poll_lineages_vec_u32(
    report: &Vec<(u32, usize)>,
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
) -> (u32, usize, usize) {
    let mut max_score = 0;
    let mut max_taxon: u32 = 0;
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report {
        let mut score = 0;
        if key.0 == (0 as u32) || key.0 == unclassified {
            continue;
        } else {
            let slice_a = &super::taxonomy_u32::get_lineage_graph(key.0, lineage_graph)[..];
            let lineage1 = Set::new(slice_a).expect("could not create set1");
            for key2 in report {
                if key2.0 == (0 as u32) || key2.0 == unclassified {
                    continue;
                } else {
                    if lineage1.contains(&key2.0) {
                        score += key2.1;
                    }
                }
            }
        }
        if score > max_score {
            max_taxon = key.0;
            max_score = score;
        }
        if score == max_score {
            max_taxon =
                super::taxonomy_u32::find_lca_u32(key.0, max_taxon, lineage_graph, unclassified);
        }
    }
    (max_taxon, max_score, kmer_length)
}

//kmer poll classification plus raw count data as output
// 1. filter hits below false prob threshold
// 2. find tophits
// 3 if only one tophit, report with unique flag
pub fn kmer_poll_plus(
    report: &hashbrown::HashMap<u32, usize>,
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
    _inverse_taxonomy: &HashMap<String, u32>,
) -> (u32, usize, usize) {
    let unclassified = (taxonomy.len() + 1) as u32;
    let report_length = report.len();
    let mut tuple = (0 as u32, 0 as usize, kmer_length);
    if report_length < 2 {
        if report_length == 0 {
            return tuple;
        } else if report_length == 1 {
            for (key, value) in report {
                if (*key == unclassified) || (*key == (0 as u32)) {
                    continue;
                } else {
                    tuple = (*key, *value, kmer_length);
                }
            }
        }
    } else {
        let lineage_scores = score_lineages(&report, &taxonomy);
        if lineage_scores.is_empty() {
            return tuple;
        } else {
            let mut count_vec: Vec<_> = lineage_scores.iter().collect();
            count_vec.sort_by(|a, b| b.1.cmp(a.1));
            //significant_hits.sort_by(|a, b| b.1.cmp(a.1));
            let first_tophit = &count_vec[0].to_owned();
            let mut top_hits = Vec::new();
            for h in &count_vec {
                if h.1 == first_tophit.1 {
                    top_hits.push(h.0.to_owned())
                }
            }
            if top_hits.len() == 1 {
                tuple = (
                    first_tophit.0.to_owned(),
                    count_vec[0].1.to_owned(),
                    kmer_length,
                );
            } else {
                tuple = (
                    super::taxonomy_u32::find_lca_vector_u32_numerical(
                        &top_hits,
                        &lineage_graph,
                        unclassified,
                    ),
                    count_vec[0].1.to_owned(),
                    kmer_length,
                );
            }
        }
    }
    tuple
}

#[inline]
pub fn vec_strings_to_string(vector_in: &[String]) -> String {
    let mut comma_separated = String::new();
    for s in vector_in {
        comma_separated.push_str(&s.to_string());
        comma_separated.push_str(",");
    }
    comma_separated.pop();
    comma_separated
}

pub struct SeqRead {
    pub id: String,  //id including >
    pub seq: String, // sequence
}

impl SeqRead {
    pub fn new() -> SeqRead {
        SeqRead {
            id: String::new(),
            seq: String::new(),
        }
    }
}

pub struct SeqReadu8 {
    pub id: String,   //id including >
    pub seq: Vec<u8>, // sequence
}

impl SeqReadu8 {
    pub fn new() -> SeqReadu8 {
        SeqReadu8 {
            id: String::new(),
            seq: Vec::new(),
        }
    }
}

pub struct SeqReadstr<'a> {
    pub id: &'a str,       //id including >
    pub seq: Vec<&'a str>, // sequence
}

impl<'a> SeqReadstr<'a> {
    pub fn new() -> SeqReadstr<'a> {
        SeqReadstr {
            id: "",
            seq: Vec::new(),
        }
    }
}

#[allow(unused_assignments)]
pub fn parallel_vec_sepia(
    vec: &[(String, String)],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    k: usize,
    m: usize,
    b: u32,
) -> std::vec::Vec<(std::string::String, u32, usize, usize)> {
    let my_db: Arc<&[u32]> = Arc::new(db);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_db = my_db.clone();
            /*if r.1.len() < k {
                (r.0.to_owned(), 0 as u32, 0 as usize, 0 as usize)
            } else {*/
            let (report, observations) = //if r.1.len() == 1 {
                    /*if m == 0 {
                        kmerize_vector_skip_n_report_phf(&r.1[0], k, b, 1, &child_db, &child_phf)
                    } else {*/
                    sliding_window_minimizers_sepia_vec(&r.1, k, m, b, &child_db);
            //}
            /*} else {
                let sequence = format!("{}{}{}", &r.1[0], "N", &r.1[1]);
                /*if m == 0 {
                    kmerize_vector_skip_n_report_phf(&sequence, k, b, 1, &child_db, &child_phf)
                } else {*/
                sliding_window_minimizers_phf(&sequence, k, m, b, &child_db, &child_phf)
                //}
            };*/
            if report.is_empty() {
                (r.0.to_owned(), 0 as u32, 0 as usize, observations)
            } else {
                /*let classification = kmer_poll_plus(
                    &report,
                    observations,
                    &lineage_graph,
                    &taxonomy,
                    &inverse_taxonomy,
                );*/
                let classification =
                    poll_lineages_vec_u32(&report, observations, &lineage_graph, &taxonomy);
                //(r.0.to_owned(), 0 as u32, 0 as usize, report.len())
                (
                    r.0.to_owned(),
                    classification.0,
                    classification.1.to_owned(),
                    classification.2.to_owned(),
                )
            }
            // }
        })
        .collect();
    out_vec
}

#[allow(unused_assignments)]
pub fn parallel_vec_sepia_pre_hll(
    vec: &[(String, String, usize)],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    k: usize,
    m: usize,
    b: u32,
) -> std::vec::Vec<(
    std::string::String,
    u32,
    usize,
    usize,
    std::vec::Vec<(u32, usize)>,
    std::vec::Vec<(u32, u64)>,
    usize,
)> {
    let my_db: Arc<&[u32]> = Arc::new(db);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_db = my_db.clone();
            let (report, observations, kmer_info) =
                sliding_window_minimizers_sepia_prehll(&r.1, k, m, b, &child_db);
            if report.is_empty() {
                (
                    r.0.to_owned(),
                    0 as u32,
                    0 as usize,
                    observations,
                    report,
                    kmer_info,
                    r.2.to_owned(),//length
                )
            } else {
                let classification =
                    poll_lineages_vec_u32(&report, observations, &lineage_graph, &taxonomy);
                (
                    r.0.to_owned(),
                    classification.0,
                    classification.1.to_owned(),
                    classification.2.to_owned(),
                    report,
                    kmer_info,
                    r.2.to_owned(),//length
                )
            }
            // }
        })
        .collect();
    out_vec
}

#[allow(unused_assignments)]
pub fn stream_fasta(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize,
    value_bits: u32,
    batch_size: usize,
    prefix: &str,
    compressed_output: bool,
) {
    //let f = File::open(filenames[0]).expect("file not found");
    //let iter1 = io::BufReader::new(f).lines();
    //let mut results = Vec::new();
    let mut results: Vec<u8> = Vec::new();
    //let mut results_counts: HashMap<u32, usize> = HashMap::new();
    let (f, _format) = niffler::from_path(&filenames[0]).expect("niffler choked");
    let iter1 = io::BufReader::with_capacity(209715200, f).byte_lines();
    //let mut l = String::with_capacity(250);
    //let mut l = String::new();
    let mut vec = Vec::with_capacity(batch_size);
    let mut sub_string = String::new();
    let search_time = SystemTime::now();
    let mut count = 0;
    let mut read_count = 0;
    let mut fasta = SeqRead::new();
    for line in iter1 {
        //while iter1.read_line(&mut l).unwrap() > 0 {
        let l = line.unwrap();
        if count == 0 {
            fasta.id = l.into_string().expect("non-UTF8 encountered in header");
        } else {
            if l[0] == b'>' {
                if !sub_string.is_empty() {
                    fasta.seq = sub_string.to_owned();
                    vec.push((fasta.id, fasta.seq));
                    fasta.id = l.into_string().expect("non-UTF8 encountered in header");
                    sub_string.clear();
                }
            } else {
                //sub_string.push_str(&l);
                sub_string.push_str(&l.into_string().expect("non-UTF8 encountered in sequenc"));
            }
        }
        count += 1;
        if vec.len() % batch_size == 0 {
            let c = parallel_vec_sepia(
                &vec,
                db,
                &taxonomy,
                &lineage_graph,
                kmer_size,
                m_size,
                value_bits,
            );
            read_count += c.len();
            eprint!(" {} reads classified\r", read_count);
            //results.append(&mut c);
            for id in c {
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
            /*for id in c {
                file.write_all(
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        id.0, id.1, id.2, id.3, id.4, id.5
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }*/
            vec.clear();
        }
        //l.clear();
    }
    //fasta.seq = sub_string;
    vec.push((fasta.id, sub_string));
    let c = parallel_vec_sepia(
        &vec,
        db,
        &taxonomy,
        &lineage_graph,
        kmer_size,
        m_size,
        value_bits,
    );
    read_count += c.len();
    //results.append(&mut c);
    for id in c {
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes());
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    eprintln!("Writing results to file...");
    if compressed_output == true {
        let mut file = BufWriter::new(GzEncoder::new(
            File::create(format!("{}_reads.gz", prefix)).expect("could not create outfile!"),
            Compression::default(),
        ));
        file.write_all(&results).expect("could not write results!");
    } else {
        std::fs::write(format!("{}_reads.txt", prefix), &results)
            .expect("Can't write results to file");
    }
}

//for benchmarking search and classification algorithms
#[allow(unused_assignments)]
pub fn per_read_stream_not_parallel(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize, //if m and k are set to the same value this effectively makes it kmer based
    value_bits: u32,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8, //if set to 0 can also be used for fasta
    compressed_output: bool,
) /* -> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    //let mut results = Vec::new();
    let mut results: Vec<u8> = Vec::new();
    //let mut results = Vec::new();
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut read_count = 0;
    while let Some(record1) = reader1.next() {
        let seqrec1 = record1.expect("invalid record in forward file");
        read_count += 1;
        //let id = str::from_utf8(seqrec1.id()).unwrap().to_string();
        let (report, observations) = sliding_window_minimizers_sepia_vec(
            &str::from_utf8(&seqrec1.raw_seq()).expect("could not read sequence"),
            kmer_size,
            m_size,
            value_bits,
            &db,
        );
        if report.is_empty() {
            results.extend(
                format!(
                    "{}\tno hits\t0\t{}\n",
                    str::from_utf8(seqrec1.id()).unwrap(),
                    observations
                )
                .as_bytes(),
            );
        } else {
            let classification =
                poll_lineages_vec_u32(&report, observations, &lineage_graph, &taxonomy);
            results.extend(
                format!(
                    "{}\t{}\t{}\t{}\n",
                    str::from_utf8(seqrec1.id()).unwrap(),
                    taxonomy[&classification.0],
                    classification.1,
                    classification.2
                )
                .as_bytes(),
            );
        }
        if read_count % 50000 == 0 {
            eprint!(" {} reads classified\r", read_count);
        }
    }
    eprint!(" {} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    eprintln!("Writing results to file...");
    if compressed_output == true {
        let mut file = BufWriter::new(GzEncoder::new(
            File::create(format!("{}_reads.gz", prefix)).expect("could not create outfile!"),
            Compression::default(),
        ));
        file.write_all(&results).expect("could not write results!");
    } else {
        std::fs::write(format!("{}_reads.txt", prefix), &results)
            .expect("Can't write results to file");
    }
}

#[allow(unused_assignments)]
pub fn per_read_stream_se(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize, //if m and k are set to the same value this effectively makes it kmer based
    value_bits: u32,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8, //if set to 0 can also be used for fasta
    hll: bool,
    compressed_output: bool,
) /* -> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    //let mut results = Vec::new();
    let mut results: Vec<u8> = Vec::new();
    let mut results_counts: HashMap<u32, usize> = HashMap::new();
    let mut results_total_length: HashMap<u32, usize> = HashMap::new();
    let mut results_ratios: HashMap<u32, f64> = HashMap::new();
    let mut hll_map: HashMap<u32, HyperLogLog<u64>> = HashMap::new();
    let mut final_counts: HashMap<u32, usize> = HashMap::new(); //actually counts of minimizers
    let mut hll_map_taxon: HashMap<u32, HyperLogLog<u64>> = HashMap::new();
    let mut final_counts_taxon: HashMap<u32, usize> = HashMap::new(); //actually counts of minimizers
    let mut vec = Vec::with_capacity(batch_size);
    //let mut results = Vec::new();
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut read_count = 0;
    while let Some(record1) = reader1.next() {
        let seqrec1 = record1.expect("invalid record in forward file");
        let seqrec1_length = seqrec1.seq().len();
        read_count += 1;
        if qual_offset == 0 {
            vec.push((
                str::from_utf8(seqrec1.id()).unwrap().to_string(),
                str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
                seqrec1_length,
            ));
        } else {
            vec.push((
                str::from_utf8(seqrec1.id()).unwrap().to_string(),
                super::seq::qual_mask(
                    &str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
                    &str::from_utf8(&seqrec1.qual().unwrap())
                        .expect("error qual reverse read")
                        .to_string(),
                    qual_offset,
                ),
                seqrec1_length,
            ))
        }
        if read_count % batch_size == 0 {
            let c = parallel_vec_sepia_pre_hll(
                &vec,
                db,
                &taxonomy,
                &lineage_graph,
                kmer_size,
                m_size,
                value_bits,
            );
            eprint!("{} reads classified\r", read_count);
            for id in c {
                *results_counts.entry(id.1).or_insert(0) += 1;
                *results_total_length.entry(id.1).or_insert(0) += id.6;
                *final_counts_taxon.entry(id.1).or_insert(0) += id.5.len();
                for c in id.4 {
                    *final_counts.entry(c.0).or_insert(0) += c.1;
                }
                if hll {
                    for f in &id.5 {
                        if hll_map_taxon.contains_key(&id.1) {
                            hll_map_taxon.get_mut(&id.1).unwrap().insert(&f.1);
                        } else {
                            let mut hllp = HyperLogLog::new(0.001);
                            hllp.insert(&f.1);
                            hll_map_taxon.insert(id.1, hllp);
                        }
                        if hll_map.contains_key(&f.0) {
                            hll_map.get_mut(&f.0).unwrap().insert(&f.1);
                        } else {
                            let mut hllp = HyperLogLog::new(0.001);
                            hllp.insert(&f.1);
                            hll_map.insert(f.0, hllp);
                        }
                    }
                }
                if id.2 == 0 {
                    *results_ratios.entry(id.1).or_insert(0.0) += 0.0;
                } else {
                    *results_ratios.entry(id.1).or_insert(0.0) += id.2 as f64 / id.3 as f64;
                }
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
            vec.clear();
        }
    }
    let c = parallel_vec_sepia_pre_hll(
        &vec,
        db,
        &taxonomy,
        &lineage_graph,
        kmer_size,
        m_size,
        value_bits,
    );
    eprint!("{} reads classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} reads in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    //results.append(&mut c);
    for id in c {
        *results_counts.entry(id.1).or_insert(0) += 1;
        *results_total_length.entry(id.1).or_insert(0) += id.6;
        *final_counts_taxon.entry(id.1).or_insert(0) += id.5.len();
        for c in id.4 {
            *final_counts.entry(c.0).or_insert(0) += c.1;
        }
        if hll {
            for f in &id.5 {
                if hll_map_taxon.contains_key(&id.1) {
                    hll_map_taxon.get_mut(&id.1).unwrap().insert(&f.1);
                } else {
                    let mut hllp = HyperLogLog::new(0.001);
                    hllp.insert(&f.1);
                    hll_map_taxon.insert(id.1, hllp);
                }
                if hll_map.contains_key(&f.0) {
                    hll_map.get_mut(&f.0).unwrap().insert(&f.1);
                } else {
                    let mut hllp = HyperLogLog::new(0.001);
                    hllp.insert(&f.1);
                    hll_map.insert(f.0, hllp);
                }
            }
        }
        if id.2 == 0 {
            *results_ratios.entry(id.1).or_insert(0.0) += 0.0;
        } else {
            *results_ratios.entry(id.1).or_insert(0.0) += id.2 as f64 / id.3 as f64;
        }
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes());
    }
    if compressed_output == true {
        let mut file = BufWriter::new(GzEncoder::new(
            File::create(format!("{}_classification.gz", prefix))
                .expect("could not create outfile!"),
            Compression::default(),
        ));
        file.write_all(&results).expect("could not write results!");
    } else {
        std::fs::write(format!("{}_classification.txt", prefix), &results)
            .expect("Can't write results to file");
    }
    let mut file =
        File::create(format!("{}_summary.txt", prefix)).expect("could not create outfile!");
    for (key, value) in results_counts {
        if final_counts.contains_key(&key) {
            if hll {
                if hll_map_taxon.contains_key(&key) {
                    file.write_all(
                        format!(
                            "{}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}\t{:.3}\n",
                            taxonomy[&key],
                            value,
                            results_ratios[&key] / value as f64,
                            final_counts[&key],
                            hll_map[&key].len().trunc() as usize,
                            final_counts[&key] as f64 / hll_map[&key].len(),
                            hll_map_taxon[&key].len().trunc() as usize,
                            final_counts_taxon[&key] as f64 / hll_map_taxon[&key].len(),
                        )
                        .as_bytes(),
                    )
                    .expect("could not write results!");
                } else {
                    file.write_all(
                        format!(
                            "{}\t{}\t{:.3}\t{}\t{}\t{:.3}\tNA\tNA\n",
                            taxonomy[&key],
                            value,
                            results_ratios[&key] / value as f64,
                            final_counts[&key],
                            hll_map[&key].len().trunc() as usize,
                            final_counts[&key] as f64 / hll_map[&key].len(),
                        )
                        .as_bytes(),
                    )
                    .expect("could not write results!");
                }
            } else {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\t{}\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                        results_total_length[&key]
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
        } else {
            if hll {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\tNA\tNA\tNA\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            } else {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
        }
    }
}

//needletail version
pub fn per_read_stream_pe(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize, //0 == no m, otherwise minimizer
    value_bits: u32,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8,
    hll: bool,
    compressed_output: bool,
) {
    let search_time = SystemTime::now();
    let mut read_count = 0;
    let mut vec = Vec::with_capacity(batch_size);
    let mut results = Vec::new();
    let mut results_counts: HashMap<u32, usize> = HashMap::new();
    let mut results_ratios: HashMap<u32, f64> = HashMap::new();
    let mut results_total_length: HashMap<u32, usize> = HashMap::new();
    let mut final_counts: HashMap<u32, usize> = HashMap::new(); //actually counts of minimizers
    let mut hll_map: HashMap<u32, HyperLogLog<u64>> = HashMap::new();
    let mut hll_map_taxon: HashMap<u32, HyperLogLog<u64>> = HashMap::new();
    let mut final_counts_taxon: HashMap<u32, usize> = HashMap::new(); //actually counts of minimizers
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut reader2 = parse_fastx_file(&filenames[1]).expect("invalid path/file");
    while let Some(record1) = reader1.next() {
        let seqrec1 = record1.expect("invalid record in forward file");
        read_count += 1;
        if let Some(record2) = reader2.next() {
            let seqrec2 = record2.expect("invalid record in reverse file");
            let sec_length = seqrec1.seq().len() + seqrec2.seq().len();
            if qual_offset == 0 {
                vec.push((
                    str::from_utf8(seqrec1.id()).unwrap().to_string(),
                    str::from_utf8(&seqrec1.seq()).unwrap().to_string()
                        + &"N"
                        + str::from_utf8(&seqrec2.seq()).unwrap(),
                    sec_length,
                ));
            } else {
                vec.push((
                    str::from_utf8(seqrec1.id()).unwrap().to_string(),
                    super::seq::qual_mask(
                        &str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
                        &str::from_utf8(&seqrec1.qual().unwrap())
                            .expect("error qual reverse read")
                            .to_string(),
                        qual_offset,
                    ) + &"N"
                        + &super::seq::qual_mask(
                            &str::from_utf8(&seqrec2.seq()).unwrap().to_string(),
                            &str::from_utf8(&seqrec2.qual().unwrap())
                                .expect("error qual reverse read")
                                .to_string(),
                            qual_offset,
                        ),
                        sec_length,
                ))
            }
        }
        if vec.len() == batch_size {
            let c = parallel_vec_sepia_pre_hll(
                &vec,
                &db,
                &taxonomy,
                &lineage_graph,
                kmer_size,
                m_size,
                value_bits,
            );
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                *results_counts.entry(id.1).or_insert(0) += 1;
                *final_counts_taxon.entry(id.1).or_insert(0) += id.5.len();
                *results_total_length.entry(id.1).or_insert(0) += id.6;
                for c in id.4 {
                    *final_counts.entry(c.0).or_insert(0) += c.1;
                }
                if hll {
                    for f in &id.5 {
                        if hll_map_taxon.contains_key(&id.1) {
                            hll_map_taxon.get_mut(&id.1).unwrap().insert(&f.1);
                        } else {
                            let mut hllp = HyperLogLog::new(0.001);
                            hllp.insert(&f.1);
                            hll_map_taxon.insert(id.1, hllp);
                        }
                        if hll_map.contains_key(&f.0) {
                            hll_map.get_mut(&f.0).unwrap().insert(&f.1);
                        } else {
                            let mut hllp = HyperLogLog::new(0.001);
                            hllp.insert(&f.1);
                            hll_map.insert(f.0, hllp);
                        }
                    }
                }
                if id.2 == 0 {
                    *results_ratios.entry(id.1).or_insert(0.0) += 0.0;
                } else {
                    *results_ratios.entry(id.1).or_insert(0.0) += id.2 as f64 / id.3 as f64;
                }
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
            vec.clear();
        }
    }
    let c = parallel_vec_sepia_pre_hll(
        &vec,
        &db,
        &taxonomy,
        &lineage_graph,
        kmer_size,
        m_size,
        value_bits,
    );
    eprint!("{} read pairs classified\r", read_count);
    for id in c {
        *results_counts.entry(id.1).or_insert(0) += 1;
        *results_total_length.entry(id.1).or_insert(0) += id.6;
        *final_counts_taxon.entry(id.1).or_insert(0) += id.5.len();
        for c in id.4 {
            *final_counts.entry(c.0).or_insert(0) += c.1;
        }
        if hll {
            for f in &id.5 {
                if hll_map_taxon.contains_key(&id.1) {
                            hll_map_taxon.get_mut(&id.1).unwrap().insert(&f.1);
                        } else {
                            let mut hllp = HyperLogLog::new(0.001);
                            hllp.insert(&f.1);
                            hll_map_taxon.insert(id.1, hllp);
                        }
                if hll_map.contains_key(&f.0) {
                    hll_map.get_mut(&f.0).unwrap().insert(&f.1);
                } else {
                    let mut hllp = HyperLogLog::new(0.001);
                    hllp.insert(&f.1);
                    hll_map.insert(f.0, hllp);
                }
            }
        }
        if id.2 == 0 {
            *results_ratios.entry(id.1).or_insert(0.0) += 0.0;
        } else {
            *results_ratios.entry(id.1).or_insert(0.0) += id.2 as f64 / id.3 as f64;
        }
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes());
    }
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    if compressed_output == true {
        let mut file = BufWriter::new(GzEncoder::new(
            File::create(format!("{}_classification.gz", prefix))
                .expect("could not create outfile!"),
            Compression::default(),
        ));
        file.write_all(&results).expect("could not write results!");
    } else {
        std::fs::write(format!("{}_classification.txt", prefix), &results)
            .expect("Can't write results to file");
    }
    let mut file =
        File::create(format!("{}_summary.txt", prefix)).expect("could not create outfile!");
        for (key, value) in results_counts {
        if final_counts.contains_key(&key) {
            if hll {
                if hll_map_taxon.contains_key(&key) {
                    file.write_all(
                        format!(
                            "{}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}\t{:.3}\n",
                            taxonomy[&key],
                            value,
                            results_ratios[&key] / value as f64,
                            final_counts[&key],
                            hll_map[&key].len().trunc() as usize,
                            final_counts[&key] as f64 / hll_map[&key].len(),
                            hll_map_taxon[&key].len().trunc() as usize,
                            final_counts_taxon[&key] as f64 / hll_map_taxon[&key].len(),
                        )
                        .as_bytes(),
                    )
                    .expect("could not write results!");
                } else {
                    file.write_all(
                        format!(
                            "{}\t{}\t{:.3}\t{}\t{}\t{:.3}\tNA\tNA\n",
                            taxonomy[&key],
                            value,
                            results_ratios[&key] / value as f64,
                            final_counts[&key],
                            hll_map[&key].len().trunc() as usize,
                            final_counts[&key] as f64 / hll_map[&key].len(),
                        )
                        .as_bytes(),
                    )
                    .expect("could not write results!");
                }
            } else {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\t{}\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                        results_total_length[&key],
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
        } else {
            if hll {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\tNA\tNA\tNA\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            } else {
                file.write_all(
                    format!(
                        "{}\t{}\t{:.3}\n",
                        taxonomy[&key],
                        value,
                        results_ratios[&key] / value as f64,
                    )
                    .as_bytes(),
                )
                .expect("could not write results!");
            }
        }
    }        
}

#[allow(unused_assignments)]
pub fn per_read_stream_pe_one_file(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize, //0 == no m, otherwise minimizer
    value_bits: u32,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8,
    compressed_output: bool,
) {
    let search_time = SystemTime::now();
    let mut results: Vec<u8> = Vec::new();
    let mut vec = Vec::with_capacity(batch_size);
    let _header = "".to_string();
    let _sequence = "".to_string();
    let mut read_count = 0;
    let mut record_count = 1;
    let batch = batch_size;
    let mut fastq = super::seq::FastqU8::new();
    let mut reader = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record in forward file");
        if record_count % 2 == 1 {
            //get header and seq/qual data read F
            fastq.id = str::from_utf8(seqrec.id()).unwrap().to_string();
            fastq.seq1 = seqrec.seq().to_vec();
            fastq.qual1 = seqrec.qual().unwrap().to_vec();
        } else {
            //get seq/qual data read R, add to vec and clear temp
            //fastq.seq2 = seqrec.seq();
            //fastq.qual2 = seqrec.qual();
            read_count += 1;
            vec.push((
                fastq.id.to_owned(),
                if qual_offset == 0 {
                    let mut merged = fastq.seq1.to_owned();
                    merged.push(b'N');
                    merged.append(&mut seqrec.seq().to_vec());
                    merged
                        .into_string()
                        .expect("non-UTF8 encountered in sequence")
                } else {
                    let mut filtered =
                        super::seq::qual_mask_u8(&fastq.seq1, &fastq.qual1, qual_offset);
                    filtered.push(b'N');
                    filtered.append(&mut super::seq::qual_mask_u8(
                        &seqrec.seq().to_vec(),
                        &seqrec.qual().unwrap().to_vec(),
                        qual_offset,
                    ));
                    filtered
                        .into_string()
                        .expect("non-UTF8 encountered in sequence")
                },
            ));
        }
        record_count += 1;
        if read_count % batch == 0 {
            let c = parallel_vec_sepia(
                &vec,
                db,
                &taxonomy,
                &lineage_graph,
                kmer_size,
                m_size,
                value_bits,
            );
            vec.clear();
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
        }
    }
    let c = parallel_vec_sepia(
        &vec,
        db,
        &taxonomy,
        &lineage_graph,
        kmer_size,
        m_size,
        value_bits,
    );
    eprint!("{} read pairs classified\r", read_count);
    match search_time.elapsed() {
        Ok(elapsed) => {
            eprintln!(
                "Classified {} read pairs in {} seconds",
                read_count,
                elapsed.as_secs()
            );
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {:?}", e);
        }
    }
    for id in c {
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes());
    }
    if compressed_output == true {
        let mut file = BufWriter::new(GzEncoder::new(
            File::create(format!("{}_reads.gz", prefix)).expect("could not create outfile!"),
            Compression::default(),
        ));
        file.write_all(&results).expect("could not write results!");
    } else {
        std::fs::write(format!("{}_reads.txt", prefix), &results)
            .expect("Can't write results to file");
    }
}
