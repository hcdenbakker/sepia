use niffler;
use std::io;
use std::io::prelude::*;

pub struct Fasta {
    pub id: String,  //id including >
    pub seq: String, // sequence
}

impl Fasta {
    pub fn new() -> Fasta {
        Fasta {
            id: String::new(),
            seq: String::new(),
        }
    }
}

pub struct Fastq {
    pub id: String,    //id including >
    pub seq1: String,  // sequence 1
    pub qual1: String, // qual seq 1
    pub seq2: String,
    pub qual2: String,
}

impl Fastq {
    pub fn new() -> Fastq {
        Fastq {
            id: String::new(),
            seq1: String::new(),
            qual1: String::new(),
            seq2: String::new(),
            qual2: String::new(),
        }
    }
}

pub struct FastqU8 {
    pub id: String,     //id including >
    pub seq1: Vec<u8>,  // sequence 1
    pub qual1: Vec<u8>, // qual seq 1
    pub seq2: Vec<u8>,
    pub qual2: Vec<u8>,
}

impl FastqU8 {
    pub fn new() -> FastqU8 {
        FastqU8 {
            id: String::new(),
            seq1: Vec::new(),
            qual1: Vec::new(),
            seq2: Vec::new(),
            qual2: Vec::new(),
        }
    }
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

#[inline]
pub fn qual_mask_u8(seq: &[u8], qual: &[u8], max_quality_offset: u8) -> Vec<u8> {
    let filtered = seq.to_owned();
    if max_quality_offset == 0 {
        filtered
    } else {
        let max_quality: u8 = max_quality_offset + 33;
        /*let mut filtered =*/
        seq.iter()
            .zip(qual)
            .map(|(base, qual)| {
                if (*qual) < max_quality {
                    b'N'
                } else {
                    base.to_owned()
                }
            })
            .collect()
    }
    //filtered
}

#[inline]
pub fn qual_mask(seq: &str, qual: &str, max_quality_offset: u8) -> String {
    if max_quality_offset == 0 {
        seq.to_owned()
    } else {
        let max_quality: u8 = max_quality_offset + 33;
        seq.chars()
            .zip(qual.chars())
            .map(|(base, qual)| {
                if (qual as u8) < max_quality {
                    'N'
                } else {
                    base.to_owned()
                }
            })
            .collect()
    }
}

//from needletail https://github.com/onecodex/needletail/blob/master/src/kmer.rs
#[inline]
pub fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
    }
}
#[inline]
pub fn has_no_n(seq: &[u8]) -> bool {
    //! Determines if a sequence has any non-primary four bases
    //! characters in it
    seq.iter().all(|n| is_good_base(*n))
}

#[test]
fn can_detect_no_n() {
    assert!(has_no_n(b"AAGT"));
    assert!(!has_no_n(b"NAGT"));
}
