#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_derive;

pub mod kmer;

pub mod seq;

pub mod zeroth;

pub mod bit_magic;

pub mod taxonomy_u32;

pub mod build_index;

pub mod search_bits;

pub mod search_bits_sepia;

pub mod direct_read_write;

pub mod classify_batch;

pub mod classify_batch_sepia;
