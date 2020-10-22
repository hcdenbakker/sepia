use boomphf::*;
use hashbrown::hash_map::HashMap;

#[inline]
pub fn hashed_key(key: u32, value_bits: u32) -> u32 {
    key >> value_bits
}

#[inline]
pub fn value(value: u32, value_bits: u32) -> u32 {
    value & ((1 << value_bits) - 1)
}

//original kraken2, uses the first bits of key as shortened key
#[inline]
pub fn populate_sepia(key: u64, value: u32, value_bits: u32) -> u32 {
    let mut cell = (key >> 32) as u32;
    cell &= !((1 << value_bits) - 1);
    cell |= value;
    cell
}

//we use the last x number of bits of the original key as shortened key
#[inline]
pub fn populate(key: u64, value: u32, value_bits: u32) -> u32 {
    let mut cell = (key << value_bits) as u32;
    cell |= value;
    cell
}

#[inline]
pub fn second_hash(first_hash: u64) -> u64 {
    (first_hash >> 8) | 1
}

#[inline]
pub fn get_phf(key: &u64, table: &[u32], phf: &Mphf<u64>, value_bits: u32) -> u32 {
    //let hc = seahash::hash(&key.as_bytes());
    let x = phf.try_hash(&key);
    if x.is_none() {
        0 as u32
    } else {
        //let compacted_key = (key >> (32 + value_bits)) as u32;
        //we flip this around
        let compacted_key = ((key << value_bits) as u32) >> value_bits;
        let idx = x.unwrap() as usize;
        if value(table[idx], value_bits) == 0 {
            //we do not expect this to happen
            0 as u32
        } else if hashed_key(table[idx], value_bits) == compacted_key {
            value(table[idx], value_bits)
        } else {
            0 as u32
        }
    }
}

#[inline]
pub fn get_sepia(key: &u64, table: &[u32], capacity: usize, value_bits: u32) -> u32 {
    let hc = murmurhash3(&key);
    //let hc = seahash::hash(&key.to_ne_bytes());
    let compacted_key = (hc >> (32 + value_bits)) as u32;
    //let compacted_key = ((key << value_bits) as u32) >> value_bits;
    let mut idx = hc as usize % capacity;
    let first_idx = idx;
    let mut step = 0;
    loop {
        if value(table[idx], value_bits) == 0 {
            // value of 0 means data is 0, saves work
            break; // search over, empty cell encountered in probe
        }
        /*} else*/
        if hashed_key(table[idx], value_bits) == compacted_key {
            return value(table[idx], value_bits);
        } // else {
        if step == 0 {
            step = second_hash(hc);
        }
        idx += step as usize;
        //idx += 1; //linear probing
        idx %= capacity;
        if idx == first_idx {
            break;
        }
        //} else {
        //    break;
        // } // search over, we've exhausted the table
    }
    0 as u32
}

#[inline]
pub fn murmurhash3(key: &u64) -> u64 {
    let mut k: u64 = *key;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    k
}

pub fn compare_and_store(table: &mut Vec<u64>, key: u64, capacity: usize) -> bool {
    let hc = seahash::hash(&key.to_ne_bytes());
    let mut set_successful = false;
    let mut search_successful = false;
    let mut idx = hc as usize % capacity;
    let first_idx = idx;
    while !search_successful {
        if table[idx] == key {
            search_successful = true;
            set_successful = true;
        }
        if table[idx] == 0 as u64 {
            table[idx] = key;
            search_successful = true;
            set_successful = true;
        }
        idx += 1; //linear probing
        idx %= capacity;
        if idx == first_idx {
            break;
        }
    }
    set_successful
}

/*pub fn compare_and_set_phf_u32(
    table: &mut [u32],
    unclassified: u32,
    lineage_graph: &HashMap<u32, u32>,
    phf: &Mphf<u64>,
    finger_prints: &[u16],
    key: String, // kmer/minimizer
    new_value: u32,
    value_bits: u32,
) -> bool {
    let hc = seahash::hash(&key.as_bytes());
    let x = phf.try_hash(&hc);
    let mut set_successful = false;
    if x.is_none() {
        false
    } else {
        let idx = x.unwrap() as usize;
        if finger_prints[idx] != hc as u16 {
            return false;
        } else {
            let compacted_key = (hc >> (32 + value_bits)) as u32;
            //let unclassified = (taxonomy.len() + 1) as u32; // infer in main function!
            //let mut idx = hc as usize % capacity;
            //let idx = phf.try_hash(&hc).unwrap() as usize;
            if table[idx] == 0 as u32 {
                table[idx] = populate(hc, new_value, value_bits);
                set_successful = true;
            } else if hashed_key(table[idx], value_bits) == compacted_key {
                let old_value = value(table[idx], value_bits);
                if new_value == old_value {
                    set_successful = true
                } else if old_value == unclassified{
                        set_successful = true;
                    } else {
                        let lca = super::taxonomy_u32::find_lca_u32(
                            new_value,
                            old_value,
                            &lineage_graph,
                            unclassified,
                        );
                        table[idx] = populate(hc, lca, value_bits);
                        set_successful = true;
                    }
            } else {
                return set_successful;
                //println!("This should not happen!");
            }
        }
        set_successful
    }
}*/

pub fn compare_and_set_phf_u32(
    table: &mut [u32],
    unclassified: u32,
    lineage_graph: &HashMap<u32, u32>,
    phf: &Mphf<u64>,
    key: u64, // kmer/minimizer
    new_value: u32,
    value_bits: u32,
) -> bool {
    //let hc = seahash::hash(&key.to_ne_bytes());
    let x = phf.try_hash(&key);
    let mut set_successful = false;
    if x.is_none() {
        false
    } else {
        let idx = x.unwrap() as usize;
        //let compacted_key = (key >> (32 + value_bits)) as u32;
        //we flip this around
        let compacted_key = ((key << value_bits) as u32) >> value_bits;
        if hashed_key(table[idx], value_bits) != compacted_key {
            return false;
        } else {
            //let compacted_key = (key >> (32 + value_bits)) as u32;
            //let unclassified = (taxonomy.len() + 1) as u32; // infer in main function!
            //let mut idx = hc as usize % capacity;
            //let idx = phf.try_hash(&hc).unwrap() as usize;
            if value(table[idx], value_bits) == 0 as u32 {
                table[idx] = populate(key, new_value, value_bits);
                set_successful = true;
            } else if hashed_key(table[idx], value_bits) == compacted_key {
                let old_value = value(table[idx], value_bits);
                if new_value == old_value {
                    set_successful = true
                } else if old_value == unclassified {
                    set_successful = true;
                } else {
                    let lca = super::taxonomy_u32::find_lca_u32(
                        new_value,
                        old_value,
                        &lineage_graph,
                        unclassified,
                    );
                    table[idx] = populate(key, lca, value_bits);
                    set_successful = true;
                }
            } else {
                return set_successful;
                //println!("This should not happen!");
            }
        }
        set_successful
    }
}

pub fn compare_and_set_sepia(
    table: &mut Vec<u32>,
    unclassified: u32,
    lineage_graph: &HashMap<u32, u32>,
    key: u64,
    new_value: u32,
    value_bits: u32,
    capacity: usize,
) -> bool {
    //let hc = seahash::hash(&key.to_ne_bytes());
    let hc = murmurhash3(&key);
    let compacted_key = (hc >> (32 + value_bits)) as u32;
    //let compacted_key = ((key << value_bits) as u32) >> value_bits;
    let mut set_successful = false;
    let mut search_successful = false;
    let mut idx = hc as usize % capacity;
    let first_idx = idx;
    let mut step = 0;
    while !search_successful {
        if table[idx] == 0 as u32 {
            table[idx] = populate_sepia(hc, new_value, value_bits);
            search_successful = true;
            set_successful = true;
        }
        /*else*/
        if hashed_key(table[idx], value_bits) == compacted_key {
            let old_value = value(table[idx], value_bits);
            if new_value == old_value {
                search_successful = true;
                set_successful = true
            } else {
                if old_value == unclassified {
                    search_successful = true;
                } else {
                    let lca = super::taxonomy_u32::find_lca_u32(
                        new_value,
                        old_value,
                        &lineage_graph,
                        unclassified,
                    );
                    table[idx] = populate_sepia(hc, lca, value_bits);
                    search_successful = true;
                    set_successful = true;
                }
            }
        } //else {
        if step == 0 {
            step = second_hash(hc) as usize;
        }
        idx += step;
        idx %= capacity;
        if idx == first_idx {
            break;
        }
        //idx += 1;
        //idx %= capacity;
        //if idx == first_idx {
        //    break;
        //}
        // }
    }
    return set_successful;
}
