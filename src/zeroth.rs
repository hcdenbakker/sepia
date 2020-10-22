use super::kmer;
use rayon::prelude::*;
use std::sync::Mutex;

#[allow(unused_assignments)]
pub fn zeroth(entries: &[Vec<String>], k: usize, m: usize, c: usize) -> usize {
    let mut count = 0;
    let mini_set: std::collections::HashSet<u64> = std::collections::HashSet::default();
    let mini_set = Mutex::new(mini_set);
    let mut vec = Vec::with_capacity(c);
    for entry in entries {
        vec.push(entry[0].to_owned());
        if vec.len() % c == 0 {
            //let mut out_vec: Vec<_> = vec![];
            //out_vec = vec
            vec.par_iter_mut().for_each(|r| {
                let fasta_vec = kmer::read_fasta(r.to_owned());
                /*if m == 0 {
                    let mmers = kmer_vector(&fasta_vec, k);
                    let mut mini_set = mini_set.lock().unwrap();
                    for hash in mmers {
                        mini_set.insert(hash);
                    }
                } else {*/
                let mmers = kmer::sliding_window_numerical_zeroth(&fasta_vec, k, m);
                let mut mini_set = mini_set.lock().unwrap();
                for hash in mmers {
                    mini_set.insert(hash);
                }
                //}
            });
            //.collect();
            /*for mmers in out_vec {
                for hash in mmers {
                    mini_set.insert(hash);
                }
                count += 1;
            }*/
            count += vec.len();
            eprint!("{}/{} accessions processed.\r", count, entries.len());
            vec.clear();
        }
    }
    //let mut out_vec: Vec<_> = vec![];
    //out_vec = vec
    vec.par_iter_mut().for_each(|r| {
        let fasta_vec = kmer::read_fasta(r.to_owned());
        /*if m == 0 {
            let mmers = kmer_vector(&fasta_vec, k);
            let mut mini_set = mini_set.lock().unwrap();
            for hash in mmers {
                mini_set.insert(hash);
            }
        } else {*/
        let mmers = kmer::sliding_window_numerical_zeroth(&fasta_vec, k, m);
        let mut mini_set = mini_set.lock().unwrap();
        for hash in mmers {
            mini_set.insert(hash);
        }
        //}
    });
    /*.collect();
    for mmers in out_vec {
        for m in mmers {
            mini_set.insert(m);
        }
        count += 1;
    }*/
    count += vec.len();
    eprintln!("{} accessions processed.", count);
    vec.clear();
    let mini_set = mini_set.into_inner().unwrap();
    mini_set.len() * 256
}
