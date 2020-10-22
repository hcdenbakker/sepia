use boomphf::*;
use fasthash;
use rayon::prelude::*;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use taxon_boom::kmer;
use redis::{transaction, Commands, PipelineCommands};
use redis;
//use redis::commands::Commands;

pub fn tab_to_map(filename: String) -> std::collections::HashMap<std::string::String, Vec<String>> {
    let mut map = std::collections::HashMap::default();
    let f = File::open(filename).expect("reference file not found");
    for line in BufReader::new(f).lines() {
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

fn main() {
    println!("Hello, world!, Let's estimate!");
    let args: Vec<String> = env::args().collect();
    let map = tab_to_map(args[1].to_owned());
    //let set = HashSet::new();
    let kmer_size = args[2].parse::<usize>().unwrap();
    let m_size = args[3].parse::<usize>().unwrap();
    let mut count = 1;
    let client = redis::Client::open("redis://127.0.0.1/").expect("couldn't establish client");
    let mut con = client.get_connection().expect("couldn't connect to server");
    //let mut pipe = redis::pipe();
    //let mut mini_set: std::collections::HashSet<u64> = std::collections::HashSet::default();
    let mut vec = Vec::with_capacity(100);
    for (key, v) in map {
        vec.push(v[0].to_owned());
        if vec.len() % 100 == 0 {
            let mut out_vec: Vec<_> = vec![];
            out_vec = vec
                .par_iter()
                .map(|r| {
                    let fasta_vec = kmer::read_fasta(r.to_owned());
                    let mut kmers = fnv::FnvHashSet::default();
                    for f in fasta_vec {
                        let mmer = kmer::sliding_window_minimizers(&f, kmer_size, m_size);
                        for m in mmer {
                            kmers.insert(seahash::hash(&m.as_bytes()));
                        }
                    }
                    kmers
                })
                .collect();
            for mmers in out_vec {
                for hash in mmers {
                    let _ : () = redis::cmd("SADD").arg("my_set").arg(hash).query(& mut con).expect("ouch!");
                    //redis::cmd("SADD").arg("my_set").arg(hash);
                    //let _ : () = con.set(hash).unwrap();
                }
                eprintln!(
                    " {} taxa added",
                    count,
                );
                count += 1;
            }
            vec.clear();
        }
    }
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let fasta_vec = kmer::read_fasta(r.to_owned());
            let mut kmers = fnv::FnvHashSet::default();
            for f in fasta_vec {
                let mmer = kmer::sliding_window_minimizers(&f, kmer_size, m_size);
                for m in mmer {
                    kmers.insert(seahash::hash(&m.as_bytes()));
                }
            }
            kmers
        })
        .collect();
    for mmers in out_vec {
        for hash in mmers {
            let _ : () = redis::cmd("SADD").arg("my_set").arg(hash).query(& mut con).expect("ouch!");
            //redis::cmd("SADD").arg("my_set").arg(hash);    
            //let _ : () = con.set(hash).unwrap();
        }
        eprintln!(
            " {} taxa added",
            count,
        );
        count += 1;
    }
    vec.clear();
    eprintln!("collecting...");
    let mut iter : redis::Iter<isize> = redis::cmd("SSCAN").arg("my_set").cursor_arg(0).clone().iter(&mut con).expect("boo!"); 
    let mut keys = Vec::new(); 
    for k in iter{
        keys.push(k);
    }
    //let keys : Vec<u64> = con.smembers("my_set").expect("boo!");
    eprintln!("{} minimers", keys.len());
}
