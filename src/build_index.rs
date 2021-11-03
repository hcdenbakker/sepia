use super::bit_magic;
use super::kmer;
use super::zeroth::zeroth;
use bincode::deserialize_from;
use boomphf::*;

use fnv;
use hashbrown::HashMap;
use rayon::prelude::*;
use std;

//use std::collections::HashSet;

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::process;
use std::sync::Mutex;
use std::time::Instant;
use sysinfo::{System, SystemExt};

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct Parameters {
    pub k_size: usize,
    pub m_size: usize,
    pub value_bits: u32,
    pub lineage_graph: HashMap<u32, u32>,
    pub mode: String,
    pub taxonomy: HashMap<u32, String>
}

pub fn read_parameters_phf(path: &str) -> Parameters {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: Parameters = deserialize_from(&mut reader).expect("can't deserialize");
    deserialized
}

pub fn tab_to_map(filename: String) -> fnv::FnvHashMap<std::string::String, Vec<String>> {
    let mut map = fnv::FnvHashMap::default();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        if v.len() == 2 {
            map.insert(String::from(v[0]), vec![String::from(v[1])]);
        } else {
            map.insert(
                String::from(v[0]),
                vec![String::from(v[1]), String::from(v[2])],
            );
        }
    }
    map
}

pub fn accession_to_tax(filename: String) -> fnv::FnvHashMap<std::string::String, String> {
    let mut map = fnv::FnvHashMap::default();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        map.insert(String::from(v[0]), String::from(v[1]));
    }
    map
}

pub fn tab_to_vec(filename: String) -> Vec<Vec<String>> {
    let mut accession_vector = Vec::new();
    let f = File::open(filename).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split('\t').collect();
        let mut string_vector = Vec::new();
        for s in v {
            string_vector.push(String::from(s));
        }
        accession_vector.push(string_vector);
        /*if v.len() == 2 {
            map.insert(String::from(v[0]), vec![String::from(v[1])]);
        } else {
            map.insert(
                String::from(v[0]),
                vec![String::from(v[1]), String::from(v[2])],
            );
        }*/
    }
    accession_vector
}

#[allow(unused_assignments)]
pub fn build_vector_parallel_mutex(
    map: &[Vec<String>],
    vector_size: usize,
    k_size: usize,
    m_size: usize,
    batch: usize,
) -> Vec<u64> {
    let map_length = map.len();
    let mut processed = 0;
    //let mut taxon_vector: Vec<Vec<u16>> = vec![vec![]; 10997000000];
    //let mut hmap = HashMap::with_hasher(BuildHasherDefault<NoHashHasher>{});
    let db: Vec<u64> = vec![0; vector_size];
    let db = Mutex::new(db);
    let mut vec = Vec::with_capacity(batch);
    for accession in map {
        vec.push(accession[0].to_owned());
        if vec.len() % batch == 0 {
            /*let mut out_vec: Vec<_> = vec![];
            out_vec = vec*/
            vec.par_iter_mut().for_each(|r| {
                let fasta_vec = kmer::read_fasta(r.to_owned()); //? to_owned? make this a reference?
                for f in fasta_vec {
                    if f.len() < k_size {
                        continue;
                    } else {
                        let mmers = /*if m_size == 0 {
                            kmer::kmerize_string_set(&f, k_size)
                        } else {*/
                            kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                        //};
                        let mut db = db.lock().unwrap();
                        for m in mmers {
                            bit_magic::compare_and_store(&mut db, m, vector_size);
                        }
                        //let kmers = kmers.union(&mmers);
                    }
                }
                //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
                //kmers
            });
            //.collect();
            /*for kmers in out_vec {
                //eprintln!("Adding {} to database", a.0);
                for k in &kmers {
                    let _set = bit_magic::compare_and_store(&mut db, k.to_string(), vector_size);
                } //
            }*/

            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    //let mut out_vec: Vec<_> = vec![];
    //out_vec = vec
    vec.par_iter_mut().for_each(|r| {
        let fasta_vec = kmer::read_fasta(r.to_owned());
        for f in fasta_vec {
            if f.len() < k_size {
                continue;
            } else {
                let mmers = /*if m_size == 0 {
                    kmer::kmerize_string_set(&f, k_size)
                } else {*/
                    kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                //};
                let mut db = db.lock().unwrap();
                for m in mmers {
                    bit_magic::compare_and_store(&mut db, m, vector_size);
                }
            }
        }
        //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
        //kmers
    });
    //.collect();
    /*for kmers in out_vec {
        for k in &kmers {
            let _set = bit_magic::compare_and_store(&mut db, k.to_string(), vector_size);
        } //
    }*/
    processed += vec.len();
    eprint!("processed {}/{} accessions\r", processed, map_length);

    db.into_inner().unwrap()
}

#[allow(unused_assignments)]
pub fn build_vector_parallel(
    map: &[Vec<String>],
    vector_size: usize,
    k_size: usize,
    m_size: usize,
    batch: usize,
) -> Vec<u64> {
    let map_length = map.len();
    let mut processed = 0;
    //let mut taxon_vector: Vec<Vec<u16>> = vec![vec![]; 10997000000];
    //let mut hmap = HashMap::with_hasher(BuildHasherDefault<NoHashHasher>{});
    let mut db: Vec<u64> = vec![0; vector_size];
    let mut vec = Vec::with_capacity(batch);
    for accession in map {
        vec.push(accession[0].to_owned());
        if vec.len() % batch == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.to_owned()); //? to_owned? make this a reference?
                    let mut kmers: std::collections::HashSet<u64> =
                        std::collections::HashSet::default();
                    for f in fasta_vec {
                        if f.len() < k_size {
                            continue;
                        } else {
                            let mmers =
                                kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                            for m in mmers {
                                kmers.insert(m);
                            }
                        }
                    }
                    kmers
                })
                .collect();
            for kmers in out_vec {
                //eprintln!("Adding {} to database", a.0);
                for k in kmers {
                    let _set = bit_magic::compare_and_store(&mut db, k, vector_size);
                } //
            }

            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.to_owned());
            let kmers: std::collections::HashSet<u64> = std::collections::HashSet::default();
            for f in fasta_vec {
                if f.len() < k_size {
                    continue;
                } else {
                    let mmers = kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                    let mut kmers: std::collections::HashSet<u64> =
                        std::collections::HashSet::default();
                    for m in mmers {
                        kmers.insert(m);
                    }
                }
            }
            kmers
        })
        .collect();
    for kmers in out_vec {
        for k in kmers {
            let _set = bit_magic::compare_and_store(&mut db, k, vector_size);
        }
    }
    processed += vec.len();
    eprint!("processed {}/{} accessions\r", processed, map_length);

    db
}

//create hash and 'empty' db
pub fn hash_fingerprints(
    map: &[Vec<String>],
    kmer_size: usize,
    m_size: usize,
    batch: usize,
    value_bits: u32,
    gamma: f64, //5.0 default?
) -> (Mphf<u64>, Vec<u32>) {
    eprintln!("Estimating total number of unique kmers/minimizers...");
    let accessions_estimate = zeroth(map, kmer_size, m_size, batch);
    let vector_size = (accessions_estimate as f64 / 0.9) as usize;
    eprintln!("Estimated number of minimizers: {}\nInitial size of compact hash set set to {}\nInferring unique minimizers...", accessions_estimate, vector_size);
    let mut hash_vector = build_vector_parallel(map, vector_size, kmer_size, m_size, batch);
    eprintln!("Removal of unused (0) entries...");
    hash_vector.retain(|&x| x != 0);
    eprintln!(
        "{} elements left after removal of unused (0) entries.",
        hash_vector.len()
    );
    eprintln!("Cheeking for duplicates, if found these will be removed...");
    hash_vector.sort_unstable();
    hash_vector.dedup();
    eprintln!(
        "\nNumber of minimizers for final index: {}",
        hash_vector.len()
    );
    eprintln!("Hashing...");
    let now = Instant::now();
    //let phf = Mphf::new(5.0, &hash_vector);
    let phf = Mphf::new_parallel(gamma, &hash_vector, None);
    eprintln!(
        "Hashing ended, done in {:3} seconds.",
        now.elapsed().as_millis() as f64 / 1000.0
    );
    //memory map hash_vector here to file to decrease RAM usage?
    let mut finger_prints: Vec<u32> = vec![0; hash_vector.len()];
    eprintln!(
        "Collecting fingerprints for final index: {} finger prints",
        hash_vector.len()
    );
    let mut finger_counter: usize = 0;
    let hash_vector_length = hash_vector.len();
    for h in hash_vector {
        let idx = phf.try_hash(&h).unwrap() as usize;
        finger_prints[idx] = bit_magic::populate(h, 0_u32, value_bits);
        finger_counter += 1;
        if finger_counter % 1000 == 0 {
            eprint!(
                "processed {}/{} accessions\r",
                finger_counter, hash_vector_length
            );
        }
    }
    eprintln!("\nCollecting fingerprints done.",);
    (phf, finger_prints)
}

//a mutex here seems very slow, check if no mutex makes it better
#[allow(unused_assignments)]
pub fn build_db_bits_parallel_phf_mutex(
    //map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    accessions: &[Vec<String>],
    taxonomy: &HashMap<u32, String>, // maybe change this to &str and just create it during build
    lookup: &HashMap<String, u32>,
    lineage_graph: &HashMap<u32, u32>,
    //mut db: &[u32],
    //phf: &Mphf<u64>,
    k_size: usize,
    m_size: usize,
    gamma: f64,
    batch: usize,
) -> (std::vec::Vec<u32>, boomphf::Mphf<u64>) {
    let unclassified = (taxonomy.len() - 1) as u32;
    let value_bits = (((unclassified as f64).log(10.0) / 2_f64.log(10.0)).floor() + 1.0) as u32;
    let (phf, db) = hash_fingerprints(accessions, k_size, m_size, batch, value_bits, gamma);
    let db = Mutex::new(db);
    let map_length = accessions.len();
    let mut processed = 0;
    eprintln!("Creating final index...");
    let mut vec = Vec::with_capacity(batch);
    for accession in accessions {
        let accession_lineage = lookup[&accession[1]]; //now a u32
        vec.push((accession_lineage, accession[0].to_owned()));
        if vec.len() % batch == 0 {
            //let mut out_vec: Vec<_> = vec![];
            //out_vec = vec
            vec.par_iter_mut().for_each(|r| {
                let accession_lineage = r.0;
                let fasta_vec = kmer::read_fasta(r.1.to_owned());
                for f in fasta_vec {
                    if f.len() < k_size {
                        continue;
                    } else {
                        let mmers = /*if m_size == 0 {
                            kmer::kmerize_string_set(&f, k_size)
                        } else {*/
                            kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                        //};
                        let mut db = db.lock().unwrap();
                        for m in mmers {
                            super::bit_magic::compare_and_set_phf_u32(
                                &mut db,
                                //&taxonomy,
                                //&lookup,
                                unclassified,
                                lineage_graph,
                                &phf,
                                m,
                                accession_lineage,
                                value_bits,
                            );
                        }
                        //let kmers = kmers.union(&mmers);
                    }
                }
                //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
                //(r.0.to_owned(), kmers)
            });
            /*.collect();
            for a in out_vec {
                let accession_lineage = a.0;
                let kmers = a.1;
                for k in &kmers {
                    let _set = super::bit_magic::compare_and_set_phf_u32(
                        &mut db,
                        //&taxonomy,
                        //&lookup,
                        unclassified,
                        &lineage_graph,
                        &phf,
                        k.to_string(),
                        accession_lineage,
                        value_bits,
                    );
                } //
            }*/

            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    //let mut out_vec: Vec<_> = vec![];
    //out_vec = vec
    vec.par_iter_mut().for_each(|r| {
        let fasta_vec = kmer::read_fasta(r.1.to_owned());
        let accession_lineage = r.0;
        //let mut kmers: std::collections::HashSet<std::string::String> = HashSet::default();
        for f in fasta_vec {
            if f.len() < k_size {
                continue;
            } else {
                let mmers = /*if m_size == 0 {
                    kmer::kmerize_string_set(&f, k_size)
                } else {*/
                    kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                //};
                let mut db = db.lock().unwrap();
                for m in mmers {
                    super::bit_magic::compare_and_set_phf_u32(
                        &mut db,
                        //&taxonomy,
                        //&lookup,
                        unclassified,
                        lineage_graph,
                        &phf,
                        m,
                        accession_lineage,
                        value_bits,
                    );
                }
            }
        }
        //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
        //(r.0.to_owned(), kmers)
    });
    //.collect();
    /*for a in out_vec {
        let accession_lineage = a.0;
        let kmers = a.1;
        for k in &kmers {
            let _set = super::bit_magic::compare_and_set_phf_u32(
                &mut db,
                //&taxonomy,
                //&lookup,
                unclassified,
                &lineage_graph,
                &phf,
                k.to_string(),
                accession_lineage,
                value_bits,
            );
        } //
    }*/
    processed += vec.len();
    eprintln!("processed {}/{} accessions", processed, map_length);
    let db = db.into_inner().unwrap();
    (db, phf)
}

#[allow(unused_assignments)]
pub fn build_db_bits_parallel_phf(
    //map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    accessions: &[Vec<String>],
    taxonomy: &HashMap<u32, String>, // maybe change this to &str and just create it during build
    lookup: &HashMap<String, u32>,
    lineage_graph: &HashMap<u32, u32>,
    //mut db: &[u32],
    //phf: &Mphf<u64>,
    k_size: usize,
    m_size: usize,
    gamma: f64,
    batch: usize,
) -> (std::vec::Vec<u32>, boomphf::Mphf<u64>) {
    let unclassified = (taxonomy.len() + 1) as u32;
    let value_bits = (((unclassified as f64).log(10.0) / 2_f64.log(10.0)).floor() + 1.0) as u32;
    let (phf, mut db) = hash_fingerprints(accessions, k_size, m_size, batch, value_bits, gamma);
    let map_length = accessions.len();
    let mut processed = 0;
    eprintln!("Creating final index...");
    let mut vec = Vec::with_capacity(batch);
    for accession in accessions {
        //let accession_lineage = lookup[&accession[1]]; //now a u32
        let accession_plus_root = "root;".to_owned() + &accession[1];
        let accession_lineage = lookup[&accession_plus_root]; //now a u32, original when we have complete taxonomy_string
        vec.push((accession_lineage, accession[0].to_owned()));
        if vec.len() % batch == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.1.to_owned());
                    let mut kmers: std::collections::HashSet<u64> =
                        std::collections::HashSet::default();
                    for f in fasta_vec {
                        if f.len() < k_size {
                            continue;
                        } else {
                            let mmers =
                                kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                            for m in mmers {
                                kmers.insert(m);
                            }
                        }
                    }
                    (r.0.to_owned(), kmers)
                })
                .collect();
            for a in out_vec {
                let accession_lineage = a.0;
                let kmers = a.1;
                for k in kmers {
                    let _set = super::bit_magic::compare_and_set_phf_u32(
                        &mut db,
                        //&taxonomy,
                        //&lookup,
                        unclassified,
                        lineage_graph,
                        &phf,
                        k,
                        accession_lineage,
                        value_bits,
                    );
                } //
            }

            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.1.to_owned());
            let mut kmers: std::collections::HashSet<u64> = std::collections::HashSet::default();
            for f in fasta_vec {
                if f.len() < k_size {
                    continue;
                } else {
                    let mmers = kmer::sliding_window_minimizers_skip_n_u64(&f, k_size, m_size);
                    for m in mmers {
                        kmers.insert(m);
                    }
                }
            }
            (r.0.to_owned(), kmers)
        })
        .collect();
    for a in out_vec {
        let accession_lineage = a.0;
        let kmers = a.1;
        for k in kmers {
            let _set = super::bit_magic::compare_and_set_phf_u32(
                &mut db,
                //&taxonomy,
                //&lookup,
                unclassified,
                lineage_graph,
                &phf,
                k,
                accession_lineage,
                value_bits,
            );
        } //
    }
    processed += vec.len();
    eprintln!("processed {}/{} accessions", processed, map_length);
    (db, phf)
}

#[allow(unused_assignments)]
pub fn build_db_bits_parallel_sepia(
    //map: &fnv::FnvHashMap<std::string::String, Vec<String>>,
    accessions: &[Vec<String>],
    taxonomy: &HashMap<u32, String>, // maybe change this to &str and just create it during build
    lookup: &HashMap<String, u32>,
    lineage_graph: &HashMap<u32, u32>,
    //mut db: &[u32],
    //phf: &Mphf<u64>,
    kmer_size: usize,
    m_size: usize,
    batch: usize,
) -> std::vec::Vec<u32> {
    let unclassified = (taxonomy.len() + 1) as u32;
    let value_bits = (((unclassified as f64).log(10.0) / 2_f64.log(10.0)).floor() + 1.0) as u32;
    eprintln!("Estimating total number of unique kmers/minimizers...");
    let accessions_estimate = zeroth(accessions, kmer_size, m_size, batch);
    let vector_size = (accessions_estimate as f64 / 0.7) as usize;
    println!(
        "Estimated index size: {:.3} Gb",
        (vector_size as f64 / 4.0) / 1073741824.0
    );
    let mut sys = System::new_all();
    sys.refresh_all();
    let estimated_size_index = (vector_size as f64 / 4.0) / 1048576.0; // in KB
    let ram = sys.total_memory();
    if estimated_size_index > ram as f64 {
        eprintln!(
            "RAM {} KB is not enough to build an index of {} KB; Abort!",
            ram, estimated_size_index
        );
        process::abort();
    }
    println!("Estimated number of minimizers: {}\nInitial size of compact hash set set to {}\nBuidling index...", accessions_estimate, vector_size);
    //let (phf, mut db) = hash_fingerprints(accessions, k_size, m_size, batch, value_bits);
    let mut db: Vec<u32> = vec![0; vector_size as usize];
    let map_length = accessions.len();
    let mut processed = 0;
    eprintln!("Creating final index...");
    let mut vec = Vec::with_capacity(batch);
    for accession in accessions {
        let accession_plus_root = "root;".to_owned() + &accession[1];
        let accession_lineage = lookup
            .get(&accession_plus_root)
            .unwrap_or_else(|| panic!("Could not find {} in lookup hash", &accession_plus_root));
        vec.push((accession_lineage, accession[0].to_owned()));
        if vec.len() % batch == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.1.to_owned());
                    let mut kmers: std::collections::HashSet<u64> =
                        std::collections::HashSet::default();
                    for f in fasta_vec {
                        if f.len() < kmer_size {
                            continue;
                        } else {
                            let mmers =
                                kmer::sliding_window_minimizers_skip_n_u64(&f, kmer_size, m_size);
                            for m in mmers {
                                kmers.insert(m);
                            }
                            //fast u64 mmer inference is approximately the same for sequence and
                            //its reverse complement; this emeliorates potential differences
                            //(between ~1 and ~5% of the minimizers
                            let mmers_rc = kmer::sliding_window_minimizers_skip_n_u64(
                                &kmer::revcomp(&f),
                                kmer_size,
                                m_size,
                            );
                            for m in mmers_rc {
                                kmers.insert(m);
                            }
                        }
                    }
                    (r.0.to_owned(), kmers)
                })
                .collect();
            for a in out_vec {
                let accession_lineage = a.0;
                let kmers = a.1;
                for k in kmers {
                    let _set = super::bit_magic::compare_and_set_sepia(
                        &mut db,
                        //&taxonomy,
                        //&lookup,
                        unclassified,
                        lineage_graph,
                        //&phf,
                        k,
                        accession_lineage,
                        value_bits,
                        vector_size,
                    );
                } //
            }

            processed += vec.len();
            eprint!("processed {}/{} accessions\r", processed, map_length);
            vec.clear();
        }
    }
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.1.to_owned());
            let mut kmers: std::collections::HashSet<u64> = std::collections::HashSet::default();
            for f in fasta_vec {
                if f.len() < kmer_size {
                    continue;
                } else {
                    let mmers = kmer::sliding_window_minimizers_skip_n_u64(&f, kmer_size, m_size);
                    for m in mmers {
                        kmers.insert(m);
                    }
                    let mmers_rc = kmer::sliding_window_minimizers_skip_n_u64(
                        &kmer::revcomp(&f),
                        kmer_size,
                        m_size,
                    );
                    for m in mmers_rc {
                        kmers.insert(m);
                    }
                }
            }
            (r.0.to_owned(), kmers)
        })
        .collect();
    for a in out_vec {
        let accession_lineage = a.0;
        let kmers = a.1;
        for k in kmers {
            let _set = super::bit_magic::compare_and_set_sepia(
                &mut db,
                //&taxonomy,
                //&lookup,
                unclassified,
                lineage_graph,
                //&phf,
                k,
                accession_lineage,
                value_bits,
                vector_size,
            );
        } //
    }
    processed += vec.len();
    eprintln!("processed {}/{} accessions", processed, map_length);
    db
}
