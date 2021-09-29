use hashbrown::HashMap;
use hashbrown::HashSet;
use sdset::duo::OpBuilder;
use sdset::{Set, SetBuf, SetOperation};
use std::fs::File;
use std::io;
use std::io::BufRead;

pub fn species_to_lineage(path: &str) -> HashMap<String, String> {
    let mut map = HashMap::default();
    let f = File::open(path).expect("reference file not found");
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split(';').collect();
        map.insert(String::from(v[6]), l);
    }
    map
}

pub fn taxonomy_map_level_agnostic(
    taxonomy_vector: &[String],
) -> (
    HashMap<u32, String>, // key: u32 value: taxonomy_string
    HashMap<String, u32>, // key: taxonomy_string, value: u32
    HashMap<u32, u32>,    //unidirectional graph connecting lower to higher order taxa
    Vec<(String,String)>, //vector with potentially problematic lineages 
) {
    let mut map = HashMap::default();
    let mut map_inv: HashMap<String, u32> = HashMap::default();
    let mut taxonomy_graph: HashMap<u32, u32> = HashMap::default(); //unidirectional graph connecting lower to higher order taxa
    let mut taxa_to_u32: HashMap<String, u32> = HashMap::default();
    let root: u32 = 1;
    map.insert(root, "root".to_string());
    map_inv.insert("root".to_string(), root);
    let mut taxon_counter: u32 = 1;
    let mut lineage_u32 = Vec::new();
    let mut lineage_string = Vec::new();
    let mut taxon_lineage: Vec<(String,String)> = Vec::new();
    for lineage in taxonomy_vector {
        let v: Vec<&str> = lineage.split(';').collect();
        lineage_u32.push(1 as u32);
        lineage_string.push("root".to_string());
        for taxon in v {   //this is going to be problematic is a taxon name is found in multiple lineages! Take the whole lineage as node name
            taxon_lineage.push((taxon.to_string(), lineage_string.join(";").to_string()));
            lineage_string.push(taxon.to_string());
            if taxa_to_u32.contains_key(&lineage_string.join(";")) {
                lineage_u32.push(taxa_to_u32[&lineage_string.join(";")]);
            } else {
                taxon_counter += 1;
                lineage_u32.push(taxon_counter);
                taxa_to_u32.insert(lineage_string.join(";"), taxon_counter);
                map.insert(taxon_counter, lineage_string.join(";"));
                map_inv.insert(lineage_string.join(";"), taxon_counter);
            }
        }
        for i in (1..lineage_u32.len()).rev() {
            taxonomy_graph.insert(lineage_u32[i], lineage_u32[i - 1]);
        }
        lineage_u32.clear();
        lineage_string.clear();
    }
    let unclassified = (map.len() + 1) as u32; //maybe root?
    map.insert(0 as u32, "no hits".to_string());
    map_inv.insert("no hits".to_string(), 0 as u32);
    map.insert(unclassified, "unclassified".to_string());
    map_inv.insert("unclassified".to_string(), unclassified);
    (map, map_inv, taxonomy_graph, taxon_lineage)
}


pub fn get_lineage_graph(query: u32, graph: &HashMap<u32, u32>) -> Vec<u32> {
    let mut stop = false;
    let mut lineage = vec![query];
    while stop == false {
        let value = graph.get(&lineage[lineage.len() - 1]);
        if value == None {
            stop = true;
        } else {
            lineage.push(*value.unwrap());
        }
    }
    lineage.sort_unstable();
    lineage
}

pub fn find_lca_u32(
    taxon_1: u32,
    taxon_2: u32,
    graph: &HashMap<u32, u32>,
    unclassified: u32,
) -> u32 {
    let slice_a = &get_lineage_graph(taxon_1, graph)[..];
    let slice_b = &get_lineage_graph(taxon_2, graph)[..];
    let a = Set::new(slice_a).expect("could not create set1");
    let b = Set::new(slice_b).expect("could not create set2");
    //let a = Set::new(&[1, 2, 4, 6, 7])?;
    //let b = Set::new(&[2, 3, 4, 5, 6, 7])?;

    let op = OpBuilder::new(a, b).intersection();

    let res: SetBuf<u32> = op.into_set_buf();
    if res[..].last().is_none() {
        unclassified
    } else {
        res[..].last().unwrap().to_owned()
    }
}

pub fn b_is_a_ancestor_of_a(taxon_a: u32, taxon_b: u32, graph: &HashMap<u32, u32>) -> bool {
    let slice_a = &get_lineage_graph(taxon_a, graph)[..];
    let a = Set::new(slice_a).expect("could not create set1");
    a.contains(&taxon_b)
}

pub fn find_lca_vector_u32_numerical(
    taxa: &Vec<u32>,
    graph: &HashMap<u32, u32>,
    unclassified: u32,
) -> u32 {
    let mut first = taxa[0];
    for t in taxa {
        //let  mut new = "";
        first = find_lca_u32(first, *t, graph, unclassified);
        //first = &new;
        if first == 0 as u32 {
            break;
        }
    }
    first
}
