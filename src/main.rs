use bincode::{deserialize_from, serialize};
use boomphf::*;
use clap::{App, AppSettings, Arg, SubCommand};
use hashbrown::HashMap;
use rayon::ThreadPoolBuilder;
use sepia::direct_read_write;
use sepia::seq;
use std::alloc::System;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::process;
use std::time::SystemTime;
use fs_extra::dir::get_size;
use sysinfo::{SystemExt, RefreshKind, System as sysinfo_System};
use std::io::{BufWriter, Write};

#[macro_use]
extern crate clap;

#[global_allocator]
static GLOBAL: System = System;

fn main() {
    let matches = App::new("sepia")
            .version("0.0.2.")
            .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
            .about("perfect hash index based read classifier and more")
            .setting(AppSettings::ArgRequiredElseHelp)
            .subcommand(
                SubCommand::with_name("build")
                    .about("builds an index")
                    .version("0.1")
                    .author("Henk C. den Bakker <henkcdenbakker@gmail.com>")
                    .setting(AppSettings::ArgRequiredElseHelp)
                    .arg(
                        Arg::with_name("index")
                            .short("i")
                            .long("index")
                            .required(true)
                            .takes_value(true),
                    )
                    .arg(
                        Arg::with_name("ref_file")
                            .help("Sets the reference file to use; this is a tab delimited file with two columns; the first column contains the path to the file, the second the associated lineage (see ref_demo.txt)")
                            .required(true)
                            .short("r")
                            .takes_value(true)
                            .long("refs"),
                    )
                    .arg(
                        Arg::with_name("k-mer_size")
                            .help("Sets k-mer size")
                            .required(true)
                            .short("k")
                            .takes_value(true)
                            .long("kmer"),
                    )
                    .arg(
                        Arg::with_name("batch")
                            .help("Sets size of batch of reads to be processed in parallel (default 300)")
                            .required(false)
                            .short("c")
                            .takes_value(true)
                            .long("batch"),
                    )
                    .arg(
                        Arg::with_name("minimizer")
                            .help("minimizer size, default length minimizer is 0 (no minimizers), unless indicated otherwise")
                            .required(true)
                            .short("m")
                            .default_value("0")
                            .takes_value(true)
                            .long("minimizer"),
                    )
                    .arg(
                        Arg::with_name("threads")
                            .help("number of threads to use, if not set one thread will be used")
                            .required(false)
                            .short("p")
                            .takes_value(true)
                            .long("threads"),
                    )
                    .arg(
                        Arg::with_name("gamma")
                            .help("gamma parameter used for perfect hash function in boom mode, default value 5.0")
                            .required(false)
                            .short("g")
                            .takes_value(true)
                            .long("gamma"),
                    )
                    .arg(
                        Arg::with_name("mode")
                            .help("build kraken2-like db (sepia), or a index with a perfect hash function (boom)")
                            .required(false)
                            .short("M")
                            .default_value("sepia")
                            .takes_value(true)
                            .long("mode"),
                    ),
            )
            .subcommand(
                SubCommand::with_name("classify")
                    .about("classifies reads using an lca approach")
                    .version("0.1")
                    .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                    .setting(AppSettings::ArgRequiredElseHelp)
                    .arg(
                        Arg::with_name("index")
                            .short("i")
                            .long("index")
                            .required(true)
                            .takes_value(true)
                            .help("index to be used for search"),
                            )
                    .arg(
                        Arg::with_name("query")
                            .help("query file(-s)fastq.gz")
                            .required(true)
                            .min_values(1)
                            .short("q")
                            .takes_value(true)
                            .long("query"),
                    )
                    .arg(
                        Arg::with_name("interleaved")
                            .help("Type of sequence data; pe: paired-end fastq.gz, se: single-end fastq.gz, inter: interleaved paired-end fastq.gz, fasta: reads in fasta format")
                            .required(false)
                            .short("I")
                            .long("interleaved"),
                    )
                    .arg(
                        Arg::with_name("batch")
                            .help("Sets size of batch of reads to be processed in parallel (default 50,000)")
                            .required(false)
                            .short("c")
                            .takes_value(true)
                            .long("batch"),
                    )
                    .arg(
                        Arg::with_name("threads")
                            .help("number of threads to use, if not set the maximum available number threads will be used")
                            .required(false)
                            .short("t")
                            .takes_value(true)
                            .long("threads"),
                    )
                    .arg(
                        Arg::with_name("prefix")
                            .help("prefix for output file(-s)")
                            .required(true)
                            .short("n") //running out of options here!
                            .takes_value(true)
                            .long("prefix"),
                    )
                    .arg(
                        Arg::with_name("quality")
                            .help("kmers with nucleotides below this minimum phred score will be excluded from the analyses (default 15)")
                            .required(false)
                            .short("Q")
                            .takes_value(true)
                            .long("quality"),
                    )
                    .arg(
                        Arg::with_name("hll")
                            .help("include HLL based estimates of cardinilaty of kmers and duplicity of kmers")
                            //.required(false)
                            .short("H")
                            .takes_value(false)
                            //.default_value("false")
                            .long("hll"),
                    )
                    .arg(
                        Arg::with_name("compress_output")
                            .help("gzip compress read classification file")
                            .required(false)
                            .short("z")
                            .takes_value(true)
                            .default_value("true")
                            .long("compress_output"),
                    ),
            )
            .subcommand(
            SubCommand::with_name("batch_classify")
                .about("classifies batch of samples reads")
                .version("0.2")
                .author("Henk den Bakker. <henkcdenbakker@gmail.com>")
                .setting(AppSettings::ArgRequiredElseHelp)
                .arg(
                    Arg::with_name("index")
                        .short("i")
                        .long("index")
                        .required(true)
                        .takes_value(true)
                        .help("index to be used for search"),
                        )
                .arg(
                    Arg::with_name("query")
                        .help("tab-delimited file with samples to be classified [sample_name reads1 reads2 (optional)]")
                        .required(true)
                        .min_values(1)
                        .short("q")
                        .takes_value(true)
                        .long("query"),
                )
                .arg(
                    Arg::with_name("tag")
                        .help("tag to be include in output file names ")
                        .required(true)
                        .short("T") //
                        .takes_value(true)
                        .long("tag"),
                )
                .arg(
                    Arg::with_name("batch")
                        .help("Sets size of batch of reads to be processed in parallel (default 50,000)")
                        .required(false)
                        .short("c")
                        .takes_value(true)
                        .long("batch"),
                )
                .arg(
                    Arg::with_name("threads")
                        .help("number of threads to use, if not set the maximum available number threads will be used")
                        .required(false)
                        .short("t")
                        .takes_value(true)
                        .long("threads"),
                )
                .arg(
                    Arg::with_name("quality")
                        .help("kmers with nucleotides below this minimum phred score will be excluded from the analyses (default 15)")
                        .required(false)
                        .short("Q")
                        .takes_value(true)
                        .long("quality"),
                        )
                .arg(
                        Arg::with_name("compress_output")
                            .help("gzip compress read classification file")
                            .required(false)
                            .short("z")
                            .takes_value(true)
                            .default_value("true")
                            .long("compress_output"),
                    ),
        ).get_matches();

    if let Some(matches) = matches.subcommand_matches("build") {
        println!("Ref_file : {}", matches.value_of("ref_file").unwrap());
        println!("db : {}", matches.value_of("index").unwrap());
        println!("k-mer size: {}", matches.value_of("k-mer_size").unwrap());
        println!("minimizer size: {}", matches.value_of("minimizer").unwrap());
        println!("hash mode: {}", matches.value_of("mode").unwrap());
        let kmer = value_t!(matches, "k-mer_size", usize).unwrap_or(31);
        let mmer = value_t!(matches, "minimizer", usize).unwrap_or(0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);
        let batch = value_t!(matches, "batch", usize).unwrap_or(300);
        println!("minimizers for {} samples inferred in parallel.", batch);
        let gamma = value_t!(matches, "gamma", f64).unwrap_or(5.0);
        let mode = value_t!(matches, "mode", String).unwrap_or("sepia".to_string());
        if mmer > 31 {
            eprintln!("This version of sepia uses a 64 bit presentation of a minimizer, larger minimizers cannot be used!");
            process::abort();
        }
        //hack to work around current clap bug with default values only being &str
        let index = matches.value_of("index").unwrap();
        //is going to be fasta-file \t lineage
        let map =
            sepia::build_index::tab_to_vec(matches.value_of("ref_file").unwrap().to_string());
        //get a vec with taxonomy
        let mut taxonomy = Vec::new();
        for accession in &map {
            taxonomy.push(accession[1].to_owned());
        }
        //let (lineages, species) =
        //    taxon_boom::taxonomy_u32::taxonomy_map(matches.value_of("taxonomy").unwrap());
        let (lineages, lookup, graph, taxon_lineage) =
            sepia::taxonomy_u32::taxonomy_map_level_agnostic(&taxonomy);
        //check if some taxa have different ancestral lineages, write to file those who have
        //and leave to user to fix or accept
        fs::DirBuilder::new()
                .recursive(true)
                .create("./".to_string() + index)
                .expect("could not initiate db directory");
        let f = File::create(index.to_owned() + "/taxonomy_ambiguities.txt").expect("Unable to create file");
        let mut f = BufWriter::new(f);
        let mut taxon_lineage_set: HashMap<String, HashSet<String>> = HashMap::default();
        for l in taxon_lineage{
            if taxon_lineage_set.contains_key(&l.0){
                let mut old_set = taxon_lineage_set[&l.0].to_owned();
                old_set.insert(l.1);
                taxon_lineage_set.insert(l.0, old_set);
            }else{
            let mut new_set: HashSet<String> = HashSet::default();
            new_set.insert(l.1);
            taxon_lineage_set.insert(l.0, new_set);
            }
        }
        //now check if there are taxa with multiple different ancestral lineages
        let mut ambiguous_taxa:u64 = 0;
        for (k, v) in taxon_lineage_set{
            if v.len() > 1{
                ambiguous_taxa += 1;
                write!(f, "{}:\n", k).expect("could not write to taxonomy_ambiguities.txt file!");
                for l in v{
                    write!(f,"{}/{}\n", l, k).expect("could not write to taxonomy_ambiguities.txt file!");
                }
            }
        }
        if ambiguous_taxa > 0{
            eprintln!("Putative ambiguities in taxonomy, check the 'taxonomy_ambiguities.txt' file in the index folder if it needs fixing!")
        }
        eprintln!(
            "length taxonomy: {}, root value: {}",
            lineages.len(),
            lookup["root"]
        );
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Can't initialize ThreadPoolBuilder");
        let bits_value =
            (((lineages.len() as f64).log(10.0) / 2_f64.log(10.0)).floor() + 1.0) as u32;

        println!("{} bits reserved to store taxon information.", bits_value);

        eprintln!("Creating index...");

        if mode == "boom".to_string() {
            let (db, phf) = sepia::build_index::build_db_bits_parallel_phf(
                &map, &lineages, &lookup, &graph, kmer, mmer, gamma, batch,
            );
            println!("Saving index to file.");
            direct_read_write::do_write_u32(&("./".to_string() + index + "/db_phf.dump"), &db);

            let serialized: Vec<u8> = serialize(&phf).unwrap();
            let mut writer = File::create(&("./".to_string() + index + "/phf")).unwrap();
            writer
                .write_all(&serialized)
                .expect("problems preparing serialized perfect hash function for writing");
        } else {
            let db = sepia::build_index::build_db_bits_parallel_sepia(
                &map, &lineages, &lookup, &graph, kmer, mmer, batch,
            );
            println!("Saving index to file.");
            fs::DirBuilder::new()
                .recursive(true)
                .create("./".to_string() + index)
                .expect("could not initiate db");
            direct_read_write::do_write_u32(&("./".to_string() + index + "/db_sepia.dump"), &db);
        }
        let parameters = sepia::build_index::Parameters {
            k_size: kmer,
            m_size: mmer,
            value_bits: bits_value,
            lineage_graph: graph,
            mode: mode,
        };

        let serialized: Vec<u8> = serialize(&parameters).unwrap();
        let mut writer = File::create(&("./".to_string() + index + "/parameters")).unwrap();
        writer
            .write_all(&serialized)
            .expect("problems preparing serialized parameters for writing");

        let serialized: Vec<u8> = serialize(&lineages).unwrap();
        let mut writer = File::create(&("./".to_string() + index + "/lineages")).unwrap();
        writer
            .write_all(&serialized)
            .expect("problems preparing serialized data for writing");
    }
    if let Some(matches) = matches.subcommand_matches("classify") {
        let bigsi_time = SystemTime::now();
        let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        let index = matches.value_of("index").unwrap();
        let prefix = matches.value_of("prefix").unwrap();
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(50000);
        let interleaved = matches.is_present("interleaved");
        let compress_output = matches.value_of("compress_output").unwrap_or("true");
        let gzip_output = if compress_output == "true" {
            true
        } else {
            false
        };
        let hll = matches.is_present("hll");
        eprintln!("Value of hll is {}", hll);
        //check size index and total/available RAM
        let index_size = get_size(index).unwrap(); //size in bytes
        // we have initialize without_processes to avoid interference with the thread pool builder 
        let mut sys_info = sysinfo_System::new_with_specifics(RefreshKind::everything().without_processes());
        let ram = sys_info.total_memory();
        if index_size as f64/1024.0 > ram as f64{
            eprintln!("Memory size {} KB is not enough to load an index of {:.3} KB; Abort!", ram, index_size);
            process::abort();
         } 
        //detect format; we determine if we have fasta or fastq, or throw an errer if we cannot
        //recognize the format, or if we have a mix of fasta and fastq
        let mut formats: HashSet<String> = HashSet::new();
        for f in &fq {
            formats.insert(seq::fastq_or_fasta(f));
        }
        if formats.len() > 1
            || (formats.get(&"unknown".to_string()) == Some(&"unknown".to_string()))
        {
            eprintln!("sequence format not recognized or mixed fasta and fastq data");
            process::abort();
        }
        if ((formats.get(&"fasta".to_string()) == Some(&"fasta".to_string())) && interleaved)
            || (fq.len() > 1 && interleaved)
        {
            eprintln!("Interleaved flag does not work with multiple files or fasta-format");
            process::abort();
        }
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Can't initialize ThreadPoolBuilder");

        eprintln!("Loading index and parameters...");

        let parameters =
            sepia::build_index::read_parameters_phf(&(index.to_owned() + "/parameters"));
        if parameters.mode == "boom".to_string() {
            eprintln!("Index with perfect hash function detected!");
            let db = direct_read_write::do_read_u32(&(index.to_owned() + "/db_phf.dump"));
            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/lineages")).expect("Can't open index!"),
            );
            let colors: HashMap<u32, String> =
                deserialize_from(&mut reader).expect("can't deserialize");

            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/phf")).expect("Can't open hash function!"),
            );
            let phf: Mphf<u64> =
                deserialize_from(&mut reader).expect("can't deserialize hash function");

            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if fq.len() == 2 && (formats.get(&"fastq".to_string()) == Some(&"fastq".to_string())) {
                eprintln!("Paired-end fastq data assumed");
                sepia::search_bits::per_read_stream_pe(
                    &fq,
                    &db,
                    &colors,
                    &phf,
                    &parameters.lineage_graph,
                    parameters.k_size,
                    parameters.m_size,
                    parameters.value_bits,
                    batch,
                    prefix,
                    quality,
                    gzip_output,
                )
            };
            if fq.len() == 1 && (formats.get(&"fastq".to_string()) == Some(&"fastq".to_string())) {
                if interleaved {
                    eprintln!("Paired-end interleaved fastq data assumed");
                    sepia::search_bits::per_read_stream_pe_one_file(
                        &fq,
                        &db,
                        &colors,
                        &phf,
                        &parameters.lineage_graph,
                        parameters.k_size,
                        parameters.m_size, //0 == no m, otherwise minimizer
                        parameters.value_bits,
                        batch,
                        prefix,
                        quality, // q cutoff
                        gzip_output,
                    )
                } else {
                    eprintln!("Single-end fastq data assumed");
                    sepia::search_bits::per_read_stream_se(
                        &fq,
                        &db,
                        &colors,
                        &phf,
                        &parameters.lineage_graph,
                        parameters.k_size,
                        parameters.m_size, //0 == no m, otherwise minimizer
                        parameters.value_bits,
                        batch,
                        prefix,
                        quality, // q cutoff
                        gzip_output,
                    )
                }
            };
            if formats.get(&"fasta".to_string()) == Some(&"fasta".to_string()) {
                eprintln!("Single-end fasta data assumed");
                sepia::search_bits::per_read_stream_se(
                    &fq,
                    &db,
                    &colors,
                    &phf,
                    &parameters.lineage_graph,
                    parameters.k_size,
                    parameters.m_size, //0 == no m, otherwise minimizer
                    parameters.value_bits,
                    batch,
                    prefix,
                    0,
                    gzip_output,
                )
            };
        } else {
            eprintln!("Index with default hash function assumed!");
            let db = direct_read_write::do_read_u32(&(index.to_owned() + "/db_sepia.dump"));
            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/lineages")).expect("Can't open index!"),
            );
            let colors: HashMap<u32, String> =
                deserialize_from(&mut reader).expect("can't deserialize");

            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }
            if fq.len() == 2 && (formats.get(&"fastq".to_string()) == Some(&"fastq".to_string())) {
                eprintln!("Paired-end fastq data assumed");
                sepia::search_bits_sepia::per_read_stream_pe(
                    &fq,
                    &db,
                    &colors,
                    &parameters.lineage_graph,
                    parameters.k_size,
                    parameters.m_size,
                    parameters.value_bits,
                    batch,
                    prefix,
                    quality,
                    hll,
                    gzip_output,
                )
            };
            if fq.len() == 1 && (formats.get(&"fastq".to_string()) == Some(&"fastq".to_string())) {
                if interleaved {
                    eprintln!("Paired-end interleaved fastq data assumed");
                    sepia::search_bits_sepia::per_read_stream_pe_one_file(
                        &fq,
                        &db,
                        &colors,
                        &parameters.lineage_graph,
                        parameters.k_size,
                        parameters.m_size, //0 == no m, otherwise minimizer
                        parameters.value_bits,
                        batch,
                        prefix,
                        quality, // q cutoff
                        gzip_output,
                    )
                } else {
                    eprintln!("Single-end fastq data assumed");
                    sepia::search_bits_sepia::per_read_stream_se(
                        &fq,
                        &db,
                        &colors,
                        &parameters.lineage_graph,
                        parameters.k_size,
                        parameters.m_size, //0 == no m, otherwise minimizer
                        parameters.value_bits,
                        batch,
                        prefix,
                        quality, // q cutoff
                        hll,     //use hll
                        gzip_output,
                    )
                }
            };
            if formats.get(&"fasta".to_string()) == Some(&"fasta".to_string()) {
                eprintln!("Single-end fasta data assumed");
                sepia::search_bits_sepia::per_read_stream_se(
                    &fq,
                    &db,
                    &colors,
                    &parameters.lineage_graph,
                    parameters.k_size,
                    parameters.m_size, //0 == no m, otherwise minimizer
                    parameters.value_bits,
                    batch,
                    prefix,
                    0,
                    hll, //use hll
                    gzip_output,
                )
            };
        }
    }
    if let Some(matches) = matches.subcommand_matches("batch_classify") {
        let bigsi_time = SystemTime::now();
        let batch_samples = matches.value_of("query").unwrap();
        //let fq: Vec<_> = matches.values_of("query").unwrap().collect();
        let threads = value_t!(matches, "threads", usize).unwrap_or(0);
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("Can't initialize ThreadPoolBuilder");
        let index = matches.value_of("index").unwrap();
        let quality = value_t!(matches, "quality", u8).unwrap_or(15);
        let batch = value_t!(matches, "batch", usize).unwrap_or(50000);
        let tag = matches.value_of("tag").unwrap();
        let compress_output = matches.value_of("compress_output").unwrap_or("true");
        let gzip_output = if compress_output == "true" {
            true
        } else {
            false
        };
        eprintln!("Loading index and parameters...");
        //check size index and total/available RAM
        let index_size = get_size(index).unwrap(); //size in bytes
        // we have initialize without_processes to avoid interference with the thread pool builder
        let mut sys_info = sysinfo_System::new_with_specifics(RefreshKind::everything().without_processes());
        let ram = sys_info.total_memory();
        if index_size as f64/1024.0 > ram as f64{
            eprintln!("Memory size {} KB is not enough to load an index of {:.3} KB; Abort!", ram, index_size);
            process::abort();
         }
        let parameters =
            sepia::build_index::read_parameters_phf(&(index.to_owned() + "/parameters"));
        if parameters.mode == "boom" {
            let db = direct_read_write::do_read_u32(&(index.to_owned() + "/db_phf.dump"));
            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/lineages")).expect("Can't open index!"),
            );
            let colors: HashMap<u32, String> =
                deserialize_from(&mut reader).expect("can't deserialize");

            let parameters =
                sepia::build_index::read_parameters_phf(&(index.to_owned() + "/parameters"));

            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/phf")).expect("Can't open hash function!"),
            );
            let phf: Mphf<u64> =
                deserialize_from(&mut reader).expect("can't deserialize hash function");

            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }

            sepia::classify_batch::batch_classify(
                batch_samples,
                &db,
                &colors,
                &phf,
                &parameters.lineage_graph,
                parameters.k_size,
                parameters.m_size, //0 == no m, otherwise minimizer
                parameters.value_bits,
                batch, //batch size for multi-threading
                quality,
                tag,
                gzip_output,
            );
        } else {
            let db = direct_read_write::do_read_u32(&(index.to_owned() + "/db_sepia.dump"));
            let mut reader = BufReader::new(
                File::open(&(index.to_owned() + "/lineages")).expect("Can't open index!"),
            );
            let colors: HashMap<u32, String> =
                deserialize_from(&mut reader).expect("can't deserialize");

            let parameters =
                sepia::build_index::read_parameters_phf(&(index.to_owned() + "/parameters"));

            match bigsi_time.elapsed() {
                Ok(elapsed) => {
                    eprintln!("Index loaded in {} seconds", elapsed.as_secs());
                }
                Err(e) => {
                    // an error occurred!
                    eprintln!("Error: {:?}", e);
                }
            }

            sepia::classify_batch_sepia::batch_classify(
                batch_samples,
                &db,
                &colors,
                &parameters.lineage_graph,
                parameters.k_size,
                parameters.m_size, //0 == no m, otherwise minimizer
                parameters.value_bits,
                batch, //batch size for multi-threading
                quality,
                tag,
                gzip_output,
            );
        }
    }
}
