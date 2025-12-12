use super::bit_magic::get_phf;
use super::kmer;
use boomphf::*;
use needletail::parse_fastx_file;
use super::build_index::Parameters;

use flate2::write::GzEncoder;
use flate2::Compression;

use dashmap::DashMap;
use hashbrown::HashMap;
use rayon::prelude::*;
use smallvec::SmallVec;
use std;

use std::fs::File;

use std::io::BufWriter;
use std::io::Write;
use std::str;
use std::sync::Arc;
use std::time::SystemTime;

use bstr::ByteVec;
use sdset::Set;

const TOGGLE: u64 = 0xe37e28c4271b5a2d;

/// SIMD-accelerated batch nucleotide processing for better performance
#[inline]
fn process_sequence_batch(seq: &str, k: usize, m: usize) -> impl Iterator<Item = (u64, usize)> + '_ {
    seq.as_bytes()
        .windows(m)
        .enumerate()
        .filter_map(move |(pos, window)| {
            let mut candidate: u64 = 0;
            let mut valid = true;

            // SIMD-accelerated nucleotide conversion
            let nuc_values = kmer::nuc_to_number_simd(window);

            for &val in &nuc_values {
                if val >= 4 {
                    valid = false;
                    break;
                }
                candidate = (candidate << 2) | val;
            }

            if valid {
                Some((kmer::canonical(candidate, m), pos))
            } else {
                None
            }
        })
}

#[inline]
fn sliding_window_minimizers_phf_vec(
    seq: &str,
    k: usize,
    m: usize,
    b: u32,
    db: &[u32],
    phf: &Mphf<u64>,
) -> (Vec<(u32, usize)>, usize) {
    let report: DashMap<u32, usize> = DashMap::new();
    let mut observations = 0;
    let mut window: SmallVec<[(u64, usize); 32]> = SmallVec::new(); // SmallVec for cache-friendly small windows
    let mut window_start = 0; // Track the start of the window
    let mut taxon = 0;
    let mut current_minimizer: u64 = 0;
    let mut counter = 0;
    let mut i = 1;
    let mut candidate: u64 = 0;
    let mut mask: u64 = 1;
    mask <<= m * 2;
    mask -= 1;
    let toggle = TOGGLE & mask;
    for n in seq.bytes() {
        //for j in seq[i..i + m].bytes(){
        //let new_char = kmer::nuc_to_number(&n);
        if kmer::nuc_to_number(n) < 4 {
            candidate <<= 2;
            candidate |= kmer::nuc_to_number(n);
            counter += 1;
            if counter >= m {
                candidate &= mask;
                //candidate &= seed;
                let mut minimizer = kmer::canonical(candidate, m);
                minimizer ^= toggle;
                // Remove elements from back that are larger than current minimizer
                while window.len() > window_start && window[window.len() - 1].0 > minimizer {
                    window.pop();
                }
                window.push((minimizer, i));
                // Remove elements from front that are outside the window
                while window.len() > window_start && (window[window_start].1 as isize) < i as isize - k as isize + m as isize {
                    window_start += 1;
                }
                if i >= k {
                    if current_minimizer == window[window_start].0 ^ toggle {
                        //    continue
                        *report.entry(taxon).or_insert(0) += 1;
                    } else {
                        //continue
                        taxon = get_phf(&(window[window_start].0 ^ toggle), db, phf, b);
                        *report.entry(taxon).or_insert(0) += 1;
                        current_minimizer = window[window_start].0;
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
            window_start = 0;
        }
        // }
        //if i == length {
        //    break;
        //}
        i += 1;
    }
    let count_vec: Vec<(u32, usize)> = report.into_iter().collect();
    (count_vec, observations)
}

#[inline]
pub fn search_index_lca(map: &[String], db: &[u32], bloom_size: u64) -> HashMap<u32, usize> {
    let mut report = HashMap::default();
    for k in map {
        let position = seahash::hash(k.as_bytes()) % bloom_size;
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
        if *key == (0_u32) || *key == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[key];
            for key2 in report.keys() {
                if *key2 == (0_u32) || *key2 == unclassified {
                    continue;
                } else {
                    let lineage2 = &taxonomy[key2];
                    if lineage1 == lineage2 {
                        let count = lineage_scores.entry(*key).or_insert(0);
                        *count += &report[key];
                    } else if lineage1.contains(lineage2) {
                        let count = lineage_scores.entry(*key).or_insert(0);
                        *count += &report[key2];
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
        if *key == (0_u32) || *key == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[key];
            for key2 in report.keys() {
                if *key2 == (0_u32) || *key2 == unclassified {
                    continue;
                } else {
                    let lineage2 = &taxonomy[key2];
                    //if lineage1 == lineage2 {
                    //    score += &report[&key];
                    //} else {
                    if lineage1.contains(lineage2) {
                        score += &report[key2];
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
    report: &[(u32, usize)],
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
) -> (u32, usize, usize) {
    let mut max_score = 0;
    let mut max_taxon: u32 = 0;
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report {
        let mut score = 0;
        if key.0 == (0_u32) || key.0 == unclassified {
            continue;
        } else {
            let lineage1 = &taxonomy[&key.0];
            for key2 in report {
                if key2.0 == (0_u32) || key2.0 == unclassified {
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
    report: &[(u32, usize)],
    kmer_length: usize,
    lineage_graph: &HashMap<u32, u32>,
    taxonomy: &HashMap<u32, String>,
) -> (u32, usize, usize) {
    let mut max_score = 0;
    let mut max_taxon: u32 = 0;
    let unclassified = (taxonomy.len() + 1) as u32;
    for key in report {
        let mut score = 0;
        if key.0 == (0_u32) || key.0 == unclassified {
            continue;
        } else {
            let slice_a = &super::taxonomy_u32::get_lineage_graph(key.0, lineage_graph)[..];
            let lineage1 = Set::new(slice_a).expect("could not create set1");
            for key2 in report {
                if key2.0 == (0_u32) || key2.0 == unclassified {
                    continue;
                } else if lineage1.contains(&key2.0) {
                    score += key2.1;
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
    let mut tuple = (0_u32, 0_usize, kmer_length);
    if report_length < 2 {
        if report_length == 0 {
            return tuple;
        } else if report_length == 1 {
            for (key, value) in report {
                if (*key == unclassified) || (*key == (0_u32)) {
                    continue;
                } else {
                    tuple = (*key, *value, kmer_length);
                }
            }
        }
    } else {
        let lineage_scores = score_lineages(report, taxonomy);
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
                        lineage_graph,
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
        comma_separated.push(',');
    }
    comma_separated.pop();
    comma_separated
}

pub struct SeqRead {
    pub id: String,  //id including >
    pub seq: String, // sequence
}

impl Default for SeqRead {
    fn default() -> Self {
        Self::new()
    }
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

impl Default for SeqReadu8 {
    fn default() -> Self {
        Self::new()
    }
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

impl<'a> Default for SeqReadstr<'a> {
    fn default() -> Self {
        Self::new()
    }
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
pub fn parallel_vec_phf(
    vec: &[(String, String)],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    lineage_graph: &HashMap<u32, u32>,
    phf: &Mphf<u64>,
    k: usize,
    m: usize,
    b: u32,
) -> std::vec::Vec<(std::string::String, u32, usize, usize)> {
    let my_db: Arc<&[u32]> = Arc::new(db);
    let my_phf: Arc<&Mphf<u64>> = Arc::new(phf);
    let mut out_vec: Vec<_> = vec![];
    out_vec = vec
        .par_iter()
        .map(|r| {
            let child_db = my_db.clone();
            let child_phf = my_phf.clone();
            /*if r.1.len() < k {
                (r.0.to_owned(), 0 as u32, 0 as usize, 0 as usize)
            } else {*/
            let (report, observations) = //if r.1.len() == 1 {
                    /*if m == 0 {
                        kmerize_vector_skip_n_report_phf(&r.1[0], k, b, 1, &child_db, &child_phf)
                    } else {*/
                    sliding_window_minimizers_phf_vec(&r.1, k, m, b, &child_db, &child_phf);
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
                (r.0.to_owned(), 0_u32, 0_usize, observations)
            } else {
                /*let classification = kmer_poll_plus(
                    &report,
                    observations,
                    &lineage_graph,
                    &taxonomy,
                    &inverse_taxonomy,
                );*/
                let classification =
                    poll_lineages_vec_u32(&report, observations, lineage_graph, taxonomy);
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
pub fn per_read_stream_se(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    phf: &Mphf<u64>,
    lineage_graph: &HashMap<u32, u32>,
    kmer_size: usize,
    m_size: usize, //0 == no m, otherwise minimizer
    value_bits: u32,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8,
    compressed_output: bool,
) /* -> std::collections::HashMap<std::string::String, usize>*/
{
    let search_time = SystemTime::now();
    let mut results: Vec<u8> = Vec::new();
    let mut vec = Vec::with_capacity(batch_size);
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut read_count = 0;
    while let Some(record1) = reader1.next() {
        let seqrec1 = record1.expect("invalid record in forward file");
        read_count += 1;
        if qual_offset == 0 {
            vec.push((
                str::from_utf8(seqrec1.id()).unwrap().to_string(),
                str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
            ));
        } else {
            vec.push((
                str::from_utf8(seqrec1.id()).unwrap().to_string(),
                super::seq::qual_mask(
                    &str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
                    &str::from_utf8(seqrec1.qual().unwrap())
                        .expect("error qual reverse read")
                        .to_string(),
                    qual_offset,
                ),
            ))
        }
        if read_count % batch_size == 0 {
            let c = parallel_vec_phf(
                &vec,
                db,
                taxonomy,
                lineage_graph,
                phf,
                kmer_size,
                m_size,
                value_bits,
            );
            eprint!("{} reads classified\r", read_count);
            for id in c {
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
            vec.clear();
        }
    }
    let c = parallel_vec_phf(
        &vec,
        db,
        taxonomy,
        lineage_graph,
        phf,
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
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, taxonomy[&id.1], id.2, id.3).as_bytes());
    }
    if compressed_output {
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
pub fn per_read_stream_pe(
    filenames: &[&str],
    db: &[u32],
    phf: &Mphf<u64>,
    parameters: &Parameters,
    batch_size: usize, //batch size for multi-threading
    prefix: &str,
    qual_offset: u8,
    compressed_output: bool,
) {
    let search_time = SystemTime::now();
    let mut read_count = 0;
    let mut vec = Vec::with_capacity(batch_size);
    let mut results = Vec::new();
    let mut reader1 = parse_fastx_file(&filenames[0]).expect("invalid path/file");
    let mut reader2 = parse_fastx_file(&filenames[1]).expect("invalid path/file");
    while let Some(record1) = reader1.next() {
        let seqrec1 = record1.expect("invalid record in forward file");
        read_count += 1;
        if let Some(record2) = reader2.next() {
            let seqrec2 = record2.expect("invalid record in reverse file");
            if qual_offset == 0 {
                vec.push((
                    str::from_utf8(seqrec1.id()).unwrap().to_string(),
                    str::from_utf8(&seqrec1.seq()).unwrap().to_string()
                        + "N"
                        + str::from_utf8(&seqrec2.seq()).unwrap(),
                ));
            } else {
                vec.push((
                    str::from_utf8(seqrec1.id()).unwrap().to_string(),
                    super::seq::qual_mask(
                        &str::from_utf8(&seqrec1.seq()).unwrap().to_string(),
                        &str::from_utf8(seqrec1.qual().unwrap())
                            .expect("error qual reverse read")
                            .to_string(),
                        qual_offset,
                    ) + "N"
                        + &super::seq::qual_mask(
                            &str::from_utf8(&seqrec2.seq()).unwrap().to_string(),
                            &str::from_utf8(seqrec2.qual().unwrap())
                                .expect("error qual reverse read")
                                .to_string(),
                            qual_offset,
                        ),
                ))
            }
        }
        if vec.len() == batch_size {
            let c = parallel_vec_phf(
                &vec,
                db,
                &parameters.taxonomy,
                &parameters.lineage_graph,
                phf,
                parameters.k_size,
                parameters.m_size,
                parameters.value_bits,
            );
            eprint!("{} read pairs classified\r", read_count);
            for id in c {
                results.extend(
                    format!("{}\t{}\t{}\t{}\n", id.0, parameters.taxonomy[&id.1], id.2, id.3).as_bytes(),
                );
            }
            vec.clear();
        }
    }
    let c = parallel_vec_phf(
        &vec,
        db,
        &parameters.taxonomy,
        &parameters.lineage_graph,
        phf,
        parameters.k_size,
        parameters.m_size,
        parameters.value_bits,
    );
    eprint!("{} read pairs classified\r", read_count);
    for id in c {
        results.extend(format!("{}\t{}\t{}\t{}\n", id.0, parameters.taxonomy[&id.1], id.2, id.3).as_bytes());
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
    if compressed_output {
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
pub fn per_read_stream_pe_one_file(
    filenames: &[&str],
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    phf: &Mphf<u64>,
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
            let c = parallel_vec_phf(
                &vec,
                db,
                taxonomy,
                lineage_graph,
                phf,
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
    let c = parallel_vec_phf(
        &vec,
        db,
        taxonomy,
        lineage_graph,
        phf,
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
    if compressed_output {
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
