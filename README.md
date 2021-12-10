[![Continuous Integration](https://github.com/hcdenbakker/sepia/actions/workflows/rust.yml/badge.svg)](https://github.com/hcdenbakker/sepia/actions/workflows/rust.yml)

# Sepia
Sepia is a (taxonomic) read classifier written in Rust based on (extensions) of the kraken2 algorithms (https://github.com/DerrickWood/kraken2), hence the name Sepia (just like the kraken a cephalopod, and the rust colored
pigment derived from it's ink sac). The reason I wrote Sepia is to create a software package that can be easily adapted to novel taxonomies or updates of taxonomies (e.g. GTDB), alter-
novel algorithms or data structures, and has features I would like myself (e.g., a fast batch mode).

## Installation
### Precompiled binaries
Download one of the precompiled binaries (currently only for macOS and Linux) [here](https://github.com/hcdenbakker/sepia/releases).
### Install from source:
1. Install Rust (https://www.rust-lang.org/)
2. Clone the sepia repository:
	```
	git clone https://github.com/hcdenbakker/sepia.git
	```
3. Compile sepia with `cargo` (Rust's compiler):
   ```
    cd sepia; cargo build --release
    ```

## Running Sepia
### Building an index 

To build an index you need collection of reference files (I masked low complexity regions using NCBI's DustMasker (https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/)
 in the 'masked folder') and a tab delimited file with the path to the reference file and the taxonomy string of the reference file. In the `ref_demo.txt` file you will find an example of what this
 should look like; we have reference genomes of Salmonella and Escherichia of GTDB r95, and the associated taxonomy in string format. 

```
./masked/Salmonella_enterica_masked.fasta       d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
./masked/Salmonella_enterica_C_masked.fasta     d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica_C
./masked/Salmonella_enterica_E_masked.fasta     d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica_E
./masked/Salmonella_enterica_D_masked.fasta     d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica_D
./masked/Salmonella_bongori_masked.fasta        d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella bongori
./masked/Escherichia_flexneri_masked.fasta      d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri
./masked/Escherichia_coli_masked.fasta  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
./masked/Escherichia_dysenteriae_masked.fasta   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia dysenteriae
./masked/Escherichia_coli_D_masked.fasta        d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_D
./masked/Escherichia_albertii_masked.fasta      d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia albertii
./masked/Escherichia_coli_C_masked.fasta        d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli_C
./masked/Escherichia_marmotae_masked.fasta      d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia marmotae
./masked/Escherichia_sp000208585_masked.fasta   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia sp000208585
./masked/Escherichia_fergusonii_masked.fasta    d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia fergusonii
./masked/Escherichia_sp001660175_masked.fasta   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia sp001660175
./masked/Escherichia_sp002965065_masked.fasta   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia sp002965065
```
If we have these files in place we can use the `build` subcommand to build an index.

```
./sepia build --help

builds an index

USAGE:
    sepia build [OPTIONS] --index <index> --kmer <k-mer_size> --minimizer <minimizer> --refs <ref_file>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --batch <batch>            Sets size of batch of reads to be processed in parallel (default 300)
    -g, --gamma <gamma>            gamma parameter used for perfect hash function in boom mode, default value 5.0
    -i, --index <index>            
    -k, --kmer <k-mer_size>        Sets k-mer size
    -m, --minimizer <minimizer>    minimizer size, default length minimizer is 0 (no minimizers), unless indicated
                                   otherwise [default: 0]
    -M, --mode <mode>              build kraken2-like db (sepia), or a index with a perfect hash function (boom)
                                   [default: sepia]
    -r, --refs <ref_file>          Sets the input file to use
    -p, --threads <threads>        number of threads to use, if not set one thread will be used
```
For the demo we build an index of the Salmonella and Escherichia references using a compact hash table (the data structure for the index used in kraken2), which is the default in Sepia:

`./sepia build -r ref_demo.txt -k 31 -m 21 -p 4 -c 2 -i demo_index`

This creates an index called 'demo_index'. Alternatively we can build an index if the `-M boom` option, which will create a smaller index with a perfect hash function. While the index is smaller with a perfect hash function,
the amount of RAM needed to build it is considerably larger, this is one of the experimental elements of Sepia.

In the index directory/folder you will find a text file reporting potential taxonomic incongruities called `taxonomy_ambiguities.txt`. It is up to the user to solve these incongruities if necessary, for example:
```
s__Hypnum_recurvatum:
root;k__Viridiplantae;p__Bryophyta;c__Bryopsida;o__Hypnales;f__Hypnaceae;g__Drepanium/s__Hypnum_recurvatum
root;k__Viridiplantae;p__Bryophyta;c__Bryopsida;o__Hypnales;f__Hypnaceae;g__Hypnum/s__Hypnum_recurvatum
```
In this case the species name `s__Hypnum_recurvatum` is found in two different lineages, once in the genus `g__Drepanium`, and once (as expected) in the genus `g__Hypnum`. This clearly a situation that needs to be corrected, as it may affect the accuracy of Sepia to classify reads belonging to either `g__Drepanium` or `g__Hypnum`. Here is another example that does not need correction:
```
g__Gomphus:
root;k__Metazoa;p__Arthropoda;c__Insecta;o__Odonata;f__Gomphidae/g__Gomphus
root;k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Gomphales;f__Gomphaceae/g__Gomphus
```
In this case the genus `g__Gomphus` is found in two lineages; (i) a lineage leading up to a family of Dragonflies (`f__Gomphidae`) and (ii) a lineage leading up a family of ectomycorrhizal mushrooms (`f__Gomphaceae`). This happens quite often when taxonomies that are covered by different nomenclature codes (https://en.wikipedia.org/wiki/Nomenclature_codes) are combined in one index. Taxonomies in Sepia are encoded in such a way that this should not pose a problem.   


## Classifying reads
```
./sepia classify --help

classifies reads using an lca approach

USAGE:
    sepia classify [FLAGS] [OPTIONS] --index <index> --prefix <prefix> --query <query>

FLAGS:
    -h, --help           Prints help information
    -H, --hll            include HLL based estimates of cardinilaty of kmers and duplicity of kmers
    -I, --interleaved    Type of sequence data; pe: paired-end fastq.gz, se: single-end fastq.gz, inter: interleaved
                         paired-end fastq.gz, fasta: reads in fasta format
    -V, --version        Prints version information

OPTIONS:
    -c, --batch <batch>                        Sets size of batch of reads to be processed in parallel (default 50,000)
    -z, --compress_output <compress_output>    gzip compress read classification file [default: true]
    -i, --index <index>                        index to be used for search
    -n, --prefix <prefix>                      prefix for output file(-s)
    -Q, --quality <quality>                    kmers with nucleotides below this minimum phred score will be excluded
                                               from the analyses (default 15)
    -q, --query <query>                        query file(-s)fastq.gz
    -t, --threads <threads>                    number of threads to use, if not set the maximum available number threads
                                               will be used
```
To classify reads, we use the `classify` subcommand. Thanks to needletail, sepia will automatically detect if the files are fastq or fasta files, and if they are compressed (e.g., fastq.gz). If two files are given, Sepia will
 assume the files are paired end. If you use an interleaved paired-end file, you have to use the `-I` flag, otherwise Sepia will assume it is a single end file. The `-H` can be used to implement hyperloglog, however the currently
 implementation is slow and still experimental.
To perform read classification with the previously created index with default values:
`sepia classify -q your_data_1.fastq.gz your_data_2.fastq.gz -i demo_index -n your_output_name`

This will create two files; a file summarizing the analysis ('your_output_name_summary.txt') and a file with the individual read classifications ('your_output_name__classification.gz'). Using the following command line:
 
`sort -grk2 -t $'\t' your_output_name_summary.txt |head -3`

we can see the first 3 lines of the summary file sorted by number of reads:

```
root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter jejuni;Campylobacter jejuni subsp. jejuni	121338	0.719	36401400
root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter	14249	0.820	4274700
root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter coli	13594	0.584	4078200
```
The first column shows the complete taxonomy lineage, the second the number of reads, the third the average k-mer similarity or if minimizers are used the approximate average k-mer similarity for all reads classified as the taxon in column
one, and the fourth column shows the sum of the length of all sequences classified as the taxon in column one. The last column is particularly helpful for sequence tachnologies that produces reads of variable length (e.g., Oxford Nanopore).

If we look at the first 3 lines of the  classification file:

`zcat campy_kala_k31_classification.gz |head -3`

```
READ:1	root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter jejuni;Campylobacter jejuni subsp. jejuni	304	305
READ:2	root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter jejuni;Campylobacter jejuni subsp. jejuni	379	410
READ:3	root;cellular organisms;Bacteria;Proteobacteria;delta/epsilon subdivisions;Epsilonproteobacteria;Campylobacterales;Campylobacteraceae;Campylobacter;Campylobacter jejuni;Campylobacter jejuni subsp. jejuni	97	413                                             
```
The first column shows the fastq header, the second the classification of the read, the third the number of kmers in the read supporting this classification, and the third the total number of kmers used for the classification. Note that the total 
number of kmers will be dependent on your quality settings (the `-Q` parameter); the default is a phred score of 15, this has to be set to something like `-Q 5` for sequence error heavy technologies like Oxford nanopore.

In the scripts folder of the repository you will find a python script to generate the input for Krona plots from the classification file.

## Batch mode
If we are working with large indices (e.g., > 50Gb), loading the index into RAM usually takes longer than the actual read clasification. If you want to use Sepia on a batch of sequence data it is more time efficient to use the `batch_classify` subcommand.
The input for this subcommand is a tab delimited file pointing out the prefix to be used for the output file and the path(-s) to the input file(-s):

 ```
 my_fasta	/path/to/my/file.fasta
 my_se_fastq	/path/to/my/file.fastq.gz
 my_pe_fastq	/path/to/my/file_1.fastq.gz	/path/to/my/file_2.fastq.gz
 ```
To run sepia in batch mode we now do:

`sepia batch_classify --index demo_index --query my_tab_delimited_file.txt --tag your_suffix`

The suffix will be included in the output file, I usually use it to show in the output files which index I used. Currently the batch function cannot be used on interleaved paired end files, these will be treated as single end files.  

## Help function
Using `-h` or `--help` will get you extensive documentation on sepia and its subcommands, for example this is what you get when you type `sepia` or `sepia -h`:
 ```
sepia 0.0.1.
Henk C. den Bakker <henkcdenbakker@gmail.com>
perfect hash index based read classifyer

USAGE:
    sepia [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    batch_classify    classifies batch of samples reads
    build             builds an index
    classify          classifies reads using an lca approach
    help              Prints this message or the help of the given subcommand(s)
 ```   