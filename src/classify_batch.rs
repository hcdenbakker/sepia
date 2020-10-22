use super::build_index;
use super::search_bits;
use boomphf::*;
use hashbrown::HashMap;

#[allow(unused_assignments)]
pub fn batch_classify(
    batch_samples: &str,
    db: &[u32],
    taxonomy: &HashMap<u32, String>,
    phf: &Mphf<u64>,
    lineage_graph: &HashMap<u32, u32>,
    k: usize,
    m: usize, //0 == no m, otherwise minimizer
    value_bits: u32,
    b: usize, //batch size for multi-threading
    qual_offset: u8,
    tag: &str,
    gzip_output: bool,
) {
    //tab to map to read batch
    let batch_map = build_index::tab_to_map(batch_samples.to_string());
    //start iterating through files here
    for (accession, fq_vec) in &batch_map {
        eprintln!("Classifying {}", accession);
        let prefix = format!("{}_{}", accession, tag);
        let mut fq: Vec<_> = vec![];
        fq = fq_vec.iter().map(|r| &r[..]).collect();
        if fq[0].ends_with(".gz") {
            if fq.len() > 1 {
                search_bits::per_read_stream_pe(
                    &fq,
                    &db,
                    &taxonomy,
                    &phf,
                    &lineage_graph,
                    k,
                    m, //0 == no m, otherwise minimizer
                    value_bits,
                    b,
                    &prefix,
                    qual_offset, // q cutoff
                    gzip_output,
                )
            } else {
                search_bits::per_read_stream_se(
                    &fq,
                    &db,
                    &taxonomy,
                    &phf,
                    &lineage_graph,
                    k,
                    m, //0 == no m, otherwise minimizer
                    value_bits,
                    b,
                    &prefix,
                    qual_offset, // q cutoff
                    gzip_output,
                )
            };
        } else {
            search_bits::per_read_stream_se(
                &fq,
                &db,
                &taxonomy,
                &phf,
                &lineage_graph,
                k,
                m, //0 == no m, otherwise minimizer
                value_bits,
                b,
                &prefix,
                0, // 0 for fasta or no quality filtering
                gzip_output,
            );
        }
    }
}
