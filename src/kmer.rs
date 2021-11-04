use super::seq;
use fnv;
use std;
use std::cmp::min;
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;

pub fn read_fasta(filename: String) -> Vec<String> {
    let mut f = File::open(filename).expect("file not found");
    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect("something went wrong reading the file");
    let mut vec = Vec::new();
    let mut sub_string: String = "".to_owned();
    let mut vec_raw = Vec::new();
    for line in contents.lines() {
        let l = line.to_string();
        vec_raw.push(l);
    }
    let length_raw_vec = vec_raw.len();
    let mut count_line = 0;
    for line in vec_raw {
        count_line += 1;
        if line.starts_with('>') {
            let l = sub_string.to_string();
            if !l.is_empty() {
                vec.push(l);
            }
            sub_string.clear();
        } else if count_line == length_raw_vec {
            let l = line.to_string();
            sub_string.push_str(&l);
            let l = sub_string.to_string();
            if !l.is_empty() {
                vec.push(l);
            }
        } else {
            let l = line.to_string();
            sub_string.push_str(&l);
        }
    }
    vec
}

#[inline]
pub fn kmerize_string_set(l: &str, k: usize) -> fnv::FnvHashSet<std::string::String> {
    let l = l.to_uppercase();
    let mut set = fnv::FnvHashSet::default();
    let length_l = l.len();
    let l_r = revcomp(&l);
    if length_l < k {
        return set;
    } else {
        for i in 0..l.len() - k + 1 {
            if seq::has_no_n(l[i..i + k].as_bytes()) {
                set.insert(min(&l[i..i + k], &l_r[length_l - (i + k)..length_l - i]).to_string());
            }
        }
    }
    set
}

#[inline]
pub fn kmerize_string_set_zeroth(l: &str, k: usize) -> fnv::FnvHashSet<u64> {
    let l = l.to_uppercase();
    let mut set = fnv::FnvHashSet::default();
    let length_l = l.len();
    let l_r = revcomp(&l);
    if length_l < k {
        return set;
    } else {
        for i in 0..l.len() - k + 1 {
            if seq::has_no_n(l[i..i + k].as_bytes()) {
                let hash = seahash::hash(
                    min(&l[i..i + k], &l_r[length_l - (i + k)..length_l - i])
                        .to_string()
                        .as_bytes(),
                );
                if (hash & 1023) < 4 {
                    set.insert(hash);
                }
            }
        }
    }
    set
}

//from https://docs.rs/bio/0.30.0/src/bio/alphabets/dna.rs.html#52-54
pub fn revcomp(text: &str) -> String {
    text.chars().rev().map(|a| complement(&a)).collect()
}

lazy_static! {
    static ref COMPLEMENT: [u8; 256] = {
        let mut comp = [0; 256];
        for (v, a) in comp.iter_mut().enumerate() {
            *a = v as u8;
        }
        for (&a, &b) in b"AGCTYRWSKMDVHBN".iter().zip(b"TCGARYWSMKHBDVN".iter()) {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;  // lowercase variants
        }
        comp
    };
}

fn complement(a: &char) -> char {
    let b = *a as u8;
    COMPLEMENT[b as usize] as char
}

pub fn sliding_window_minimizers(
    seq: &str,
    k: usize,
    m: usize,
) -> fnv::FnvHashSet<std::string::String> {
    let mut window: VecDeque<(&str, usize)> = VecDeque::new();
    let mut mset = std::collections::HashSet::default();
    let r_seq = revcomp(seq);
    let length = seq.len();
    for i in 0..seq.len() - m + 1 {
        let minimizer = min(&seq[i..i + m], &r_seq[length - (i + m)..length - i]);
        while !window.is_empty() && window.back().unwrap().0 > minimizer {
            window.pop_back();
        }
        window.push_back((minimizer, i));
        while (window.front().unwrap().1 as isize) < (i as isize - k as isize + m as isize) {
            window.pop_front();
        }
        if i >= k - m && seq::has_no_n(seq[i - (k - m)..i + m].as_bytes()) {
            mset.insert(window.front().unwrap().0.to_string().to_uppercase());
        }
    }
    mset
}

pub fn sliding_window_minimizersi_zeroth(
    v: &[String],
    k: usize,
    m: usize,
) -> std::collections::HashSet<u64> {
    let mut mset = std::collections::HashSet::default();
    for seq in v {
        if seq.len() < k {
            continue;
        }
        let mut window: VecDeque<(&str, usize)> = VecDeque::new();
        let r_seq = revcomp(seq);
        let length = seq.len();
        for i in 0..seq.len() - m + 1 {
            let minimizer = min(&seq[i..i + m], &r_seq[length - (i + m)..length - i]);
            while !window.is_empty() && window.back().unwrap().0 > minimizer {
                window.pop_back();
            }
            window.push_back((minimizer, i));
            while (window.front().unwrap().1 as isize) < i as isize - k as isize + m as isize {
                window.pop_front();
            }
            if i >= k - m && seq::has_no_n(seq[i - (k - m)..i + m].as_bytes()) {
                let hash = seahash::hash(
                    window
                        .front()
                        .unwrap()
                        .0
                        .to_string()
                        .to_uppercase()
                        .as_bytes(),
                );
                if (hash & 1023) < 4 {
                    mset.insert(hash);
                }
            }
        }
    }
    mset
}

const NUC_TO_NUMBER: [u64; 256] = {
    let mut nuc_to_number = [4; 256];
    nuc_to_number[b'A' as usize] = 0x00;
    nuc_to_number[b'C' as usize] = 0x01;
    nuc_to_number[b'G' as usize] = 0x02;
    nuc_to_number[b'T' as usize] = 0x03;
    nuc_to_number[b'a' as usize] = 0x00;
    nuc_to_number[b'c' as usize] = 0x01;
    nuc_to_number[b'g' as usize] = 0x02;
    nuc_to_number[b't' as usize] = 0x03;
    nuc_to_number
};

#[inline]
pub fn nuc_to_number(nuc: u8) -> u64 {
    NUC_TO_NUMBER[nuc as usize]
}

/*
#[inline(always)]
pub fn nuc_to_number(x: u8) -> u64 {
    match x {
        b'A' => 0x00,
        b'C' => 0x01,
        b'G' => 0x02,
        b'T' => 0x03,
        _ => 0x04,
    }
}*/

//from kraken2
#[inline]
fn reverse_complement(kmer: u64, kmer_size: usize) -> u64 {
    // Reverse bits (leaving bit pairs - nucleotides - intact)
    // swap consecutive pairs
    let mut kmer = ((kmer & 0xCCCCCCCCCCCCCCCC) >> 2) | ((kmer & 0x3333333333333333) << 2);
    // swap consecutive nibbles
    kmer = ((kmer & 0xF0F0F0F0F0F0F0F0) >> 4) | ((kmer & 0x0F0F0F0F0F0F0F0F) << 4);
    // swap consecutive bytes
    kmer = ((kmer & 0xFF00FF00FF00FF00) >> 8) | ((kmer & 0x00FF00FF00FF00FF) << 8);
    // swap consecutive byte pairs
    kmer = ((kmer & 0xFFFF0000FFFF0000) >> 16) | ((kmer & 0x0000FFFF0000FFFF) << 16);
    // swap halves of 64-bit word
    kmer = (kmer >> 32) | (kmer << 32);
    // Then complement
    ((!kmer) >> (8 * 8 - kmer_size * 2)) & (((1_u64) << (kmer_size * 2)) - 1)
}

//from kraken2
#[inline]
pub fn canonical(kmer: u64, kmer_size: usize) -> u64 {
    std::cmp::min(kmer, reverse_complement(kmer, kmer_size))
}

const TOGGLE: u64 = 0xe37e28c4271b5a2d;

#[inline]
pub fn sliding_window_minimizers_skip_n_u64(seq: &str, k: usize, m: usize) -> Vec<u64> {
    let mut vec: Vec<u64> = Vec::new();
    //let mut window: VecDeque<(ArrayString<[_; m]>, usize)> = VecDeque::new(); //position minimizer
    let mut window: VecDeque<(u64, usize)> = VecDeque::new(); //position minimizer
    let mut counter = 0;
    let mut i = 1;
    //let mut j = 1; //total seq counter
    let mut candidate: u64 = 0;
    let mut mask: u64 = 1;
    mask <<= m * 2;
    mask -= 1;
    let toggle = TOGGLE & mask;
    for n in seq.bytes() {
        //for j in seq[i..i + m].bytes(){
        let new_char = nuc_to_number(n);
        if new_char < 4 {
            candidate <<= 2;
            candidate |= new_char;
            counter += 1;
            if counter >= m {
                candidate &= mask;
                let mut minimizer = canonical(candidate, m);
                minimizer ^= toggle;
                while !window.is_empty() && window.back().unwrap().0 > minimizer {
                    window.pop_back(); // we pop the last one
                }
                window.push_back((minimizer, i)); // and make add a pair with the new value at the end
                while (window.front().unwrap().1 as isize) < i as isize - k as isize + m as isize {
                    window.pop_front(); // pop the first one
                }
                if i >= k {
                    //if has_no_n(&seq[i - (k - m)..i + m].as_bytes()) {
                    //we do not want to include minimers from kmers with an N
                    //change minimizer into u64
                    vec.push(window.front().unwrap().0 ^ toggle);
                    //}
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
        i += 1;
        /*j += 1;
        if j == length {
            break;
        }*/
    }
    vec
}

#[inline]
pub fn sliding_window_numerical_zeroth(
    v: &[String],
    k: usize,
    m: usize,
) -> std::collections::HashSet<u64> {
    let mut mset = std::collections::HashSet::default();
    for sequence in v {
        if sequence.len() < k {
            continue;
        }
        let seq = &sequence.to_uppercase();
        let both_strands = [seq, &revcomp(seq)];
        let mut window: VecDeque<(u64, usize)> = VecDeque::new();
        let length = seq.len();
        let mut counter = 0;
        let mut i = 1;
        let mut candidate: u64 = 0;
        let mut mask: u64 = 1;
        mask <<= m * 2;
        mask -= 1;
        let toggle = TOGGLE & mask;
        for s in &both_strands {
            for n in s.bytes() {
                //for j in seq[i..i + m].bytes(){
                let new_char = nuc_to_number(n);
                if new_char < 4 {
                    candidate <<= 2;
                    candidate |= new_char;
                    counter += 1;
                    if counter >= m {
                        candidate &= mask;
                        let mut minimizer = canonical(candidate, m);
                        minimizer ^= toggle;
                        while !window.is_empty() && window.back().unwrap().0 > minimizer {
                            window.pop_back(); // we pop the last one
                        }
                        window.push_back((minimizer, i)); // and make add a pair with the new value at the end
                        while (window.front().unwrap().1 as isize)
                            < i as isize - k as isize + m as isize
                        {
                            window.pop_front(); // pop the first one
                        }
                        if i >= k {
                            //if has_no_n(&seq[i - (k - m)..i + m].as_bytes()) {
                            //we do not want to include minimers from kmers with an N
                            //change minimizer into u64
                            //vec.push(window.front().unwrap().0);
                            let hash =
                                seahash::hash(&(window.front().unwrap().0 ^ toggle).to_ne_bytes());
                            if (hash & 1023) < 4 {
                                mset.insert(hash);
                            }
                            //}
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
                if i == length {
                    break;
                }
                i += 1;
            }
        }
    }
    mset
}
