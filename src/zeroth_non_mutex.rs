use super::kmer;

use rayon::prelude::*;

fn kmer_vector(v: &[String], k: usize) -> std::collections::HashSet<u64> {
    let mut set = std::collections::HashSet::default();
    for string in v {
        for m in kmer::kmerize_string_set_zeroth(&string, k) {
            set.insert(m);
        }
    }
    set
}

#[allow(unused_assignments)]
pub fn zeroth(entries: &[Vec<String>], k: usize, m: usize, c: usize) -> usize {
    let mut count = 1;
    let mut mini_set: std::collections::HashSet<u64> = std::collections::HashSet::default();
    let mut vec = Vec::with_capacity(c);
    for entry in entries {
        vec.push(entry[0].to_owned());
        if vec.len() % c == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.to_owned());
                    if m == 0 {
                        kmer_vector(&fasta_vec, k)
                    } else {
                        kmer::sliding_window_minimizersi_zeroth(&fasta_vec, k, m)
                    }
                })
                .collect();
            for mmers in out_vec {
                for hash in mmers {
                    mini_set.insert(hash);
                }
                count += 1;
            }
            eprint!("{} taxa processed.\r", count - 1);
            vec.clear();
        }
    }
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.to_owned());
            if m == 0 {
                kmer_vector(&fasta_vec, k)
            } else {
                kmer::sliding_window_minimizersi_zeroth(&fasta_vec, k, m)
            }
        })
        .collect();
    for mmers in out_vec {
        for m in mmers {
            mini_set.insert(m);
        }
        count += 1;
    }
    eprintln!("{} taxa processed.", count - 1);
    vec.clear();
    mini_set.len() * 256
}
