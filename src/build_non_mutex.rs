use super::bit_magic;
use super::kmer;
use super::zeroth::zeroth;
use bincode::deserialize_from;
use boomphf::*;

use fnv;
use hashbrown::HashMap;
use rayon::prelude::*;
use std;

use std::collections::HashSet;

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;

use std::time::Instant;

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct ParametersPhf {
    pub k_size: usize,
    pub m_size: usize,
    pub value_bits: u32,
    pub lineage_graph: HashMap<u32, u32>,
}

pub fn read_parameters_phf(path: &str) -> ParametersPhf {
    let mut reader = BufReader::new(File::open(path).expect("Can't open index!"));
    let deserialized: ParametersPhf = deserialize_from(&mut reader).expect("can't deserialize");
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
                    let mut kmers: std::collections::HashSet<std::string::String> =
                        HashSet::default();
                    for f in fasta_vec {
                        if f.len() < k_size {
                            continue;
                        } else {
                            let mmers = if m_size == 0 {
                                kmer::kmerize_string_set(&f, k_size)
                            } else {
                                kmer::sliding_window_minimizers(&f, k_size, m_size)
                            };
                            for m in mmers {
                                kmers.insert(m);
                            }
                            //let kmers = kmers.union(&mmers);
                        }
                    }
                    //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
                    kmers
                })
                .collect();
            for kmers in out_vec {
                //eprintln!("Adding {} to database", a.0);
                for k in &kmers {
                    let _set = bit_magic::compare_and_store(&mut db, k.to_string(), vector_size);
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
            let mut kmers: std::collections::HashSet<std::string::String> = HashSet::default();
            for f in fasta_vec {
                if f.len() < k_size {
                    continue;
                } else {
                    let mmers = if m_size == 0 {
                        kmer::kmerize_string_set(&f, k_size)
                    } else {
                        kmer::sliding_window_minimizers(&f, k_size, m_size)
                    };
                    for m in mmers {
                        kmers.insert(m);
                    }
                }
            }
            //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
            kmers
        })
        .collect();
    for kmers in out_vec {
        for k in &kmers {
            let _set = bit_magic::compare_and_store(&mut db, k.to_string(), vector_size);
        } //
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
) -> (Mphf<u64>, Vec<u32>) {
    let accessions_estimate = zeroth(&map, kmer_size, m_size, batch);
    let vector_size = (accessions_estimate as f64 / 0.9) as usize;
    eprintln!("Estimated number of minimizers: {}\nInitial size of compact hash set set to {}\nInferring unique minimizers...", accessions_estimate, vector_size);
    let mut hash_vector = build_vector_parallel(&map, vector_size, kmer_size, m_size, batch);
    hash_vector.retain(|&x| x != 0);
    hash_vector.sort_unstable();
    hash_vector.dedup();
    eprintln!(
        "\nNumber of minimizers for final index: {}",
        hash_vector.len()
    );
    eprintln!("Hashing...");
    let now = Instant::now();
    //let phf = Mphf::new(5.0, &hash_vector);
    let phf = Mphf::new_parallel(5.0, &hash_vector, None);
    eprintln!(
        "Hashing ended, done in {:3} seconds.",
        now.elapsed().as_millis() as f64 / 1000.0
    );
    let mut finger_prints: Vec<u32> = vec![0; hash_vector.len()];
    for h in hash_vector {
        let idx = phf.try_hash(&h).unwrap() as usize;
        finger_prints[idx] = bit_magic::populate(h, 0 as u32, value_bits);
    }
    (phf, finger_prints)
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
    batch: usize,
) -> (std::vec::Vec<u32>, boomphf::Mphf<u64>) {
    let unclassified = (taxonomy.len() - 1) as u32;
    let value_bits = (((unclassified as f64).log(10.0) / 2_f64.log(10.0)).floor() + 1.0) as u32;
    let (phf, mut db) = hash_fingerprints(accessions, k_size, m_size, batch, value_bits);
    let map_length = accessions.len();
    let mut processed = 0;
    eprintln!("Creating final index...");
    let mut vec = Vec::with_capacity(batch);
    for accession in accessions {
        let accession_lineage = lookup[&accession[1]]; //now a u32
        vec.push((accession_lineage, accession[0].to_owned()));
        if vec.len() % batch == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.1.to_owned());
                    let mut kmers: std::collections::HashSet<std::string::String> =
                        HashSet::default();
                    for f in fasta_vec {
                        if f.len() < k_size {
                            continue;
                        } else {
                            let mmers = if m_size == 0 {
                                kmer::kmerize_string_set(&f, k_size)
                            } else {
                                kmer::sliding_window_minimizers(&f, k_size, m_size)
                            };
                            for m in mmers {
                                kmers.insert(m);
                            }
                            //let kmers = kmers.union(&mmers);
                        }
                    }
                    //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
                    (r.0.to_owned(), kmers)
                })
                .collect();
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
            let mut kmers: std::collections::HashSet<std::string::String> = HashSet::default();
            for f in fasta_vec {
                if f.len() < k_size {
                    continue;
                } else {
                    let mmers = if m_size == 0 {
                        kmer::kmerize_string_set(&f, k_size)
                    } else {
                        kmer::sliding_window_minimizers(&f, k_size, m_size)
                    };
                    for m in mmers {
                        kmers.insert(m);
                    }
                }
            }
            //let kmers = kmer::minimerize_vector_skip_n_set(&fasta_vec, k_size, m_size, d_size);
            (r.0.to_owned(), kmers)
        })
        .collect();
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
    }
    processed += vec.len();
    eprintln!("processed {}/{} accessions", processed, map_length);

    (db, phf)
}
