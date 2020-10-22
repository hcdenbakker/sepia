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

pub fn taxonomy_map(path: &str) -> (HashMap<u32, String>, HashMap<String, u32>) {
    let mut map = HashMap::default();
    let mut species_lineage = HashMap::default();
    let mut set: HashSet<String> = HashSet::default();
    //get sets for every taxonomic level
    let f = File::open(path).expect("reference file not found");
    let _counter: u16 = 1; //start with 1, reserve 0 for no hits/data
    for line in io::BufReader::new(f).lines() {
        let l = line.unwrap();
        let v: Vec<&str> = l.split(';').collect();
        set.insert(v.join(";"));
        set.insert(v[..6].join(";"));
        set.insert(v[..5].join(";"));
        set.insert(v[..4].join(";"));
        set.insert(v[..3].join(";"));
        set.insert(v[..2].join(";"));
        set.insert(v[0].to_string());
    }
    let mut counter: u32 = 1; //start with 1, reserve 0 for no hits/data
    map.insert(0 as u32, "no hits".to_string()); // insert zero for no hits
    for t in &set {
        let v: Vec<&str> = t.split(';').collect();
        if v.len() == 7 {
            species_lineage.insert(v[6].replace(" ", "_"), counter);
        }
        map.insert(counter, t.to_owned());
        counter += 1;
    }
    (map, species_lineage)
}

pub fn taxonomy_map_level_agnostic(
    taxonomy_vector: &[String],
) -> (
    HashMap<u32, String>, // k: u32 v: taxonomy_string
    HashMap<String, u32>, // k: taxonomy_string, v: u32
    HashMap<u32, u32>,    //unidirectional graph connecting lower to higher order taxa
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
    for lineage in taxonomy_vector {
        let v: Vec<&str> = lineage.split(';').collect();
        lineage_u32.push(1 as u32);
        lineage_string.push("root".to_string());
        for taxon in v {
            if taxa_to_u32.contains_key(taxon) {
                lineage_u32.push(taxa_to_u32[taxon]);
                lineage_string.push(taxon.to_string());
            } else {
                taxon_counter += 1;
                lineage_u32.push(taxon_counter);
                lineage_string.push(taxon.to_string());
                taxa_to_u32.insert(taxon.to_string(), taxon_counter);
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
    (map, map_inv, taxonomy_graph)
}

//to include 9 levels (subspecies and strain potentially); number taxa so higher
// order taxa always have a u32 < than compartively lower order taxa
pub fn taxonomy_map_plus(
    taxonomy_vector: &[String],
) -> (
    HashMap<u32, String>,
    HashMap<String, u32>,
    HashMap<u32, u32>,
) {
    //set per level

    let mut kingdom: HashSet<String> = HashSet::default();
    let mut phylum: HashSet<String> = HashSet::default();
    let mut class: HashSet<String> = HashSet::default();
    let mut order: HashSet<String> = HashSet::default();
    let mut family: HashSet<String> = HashSet::default();
    let mut genus: HashSet<String> = HashSet::default();
    let mut species: HashSet<String> = HashSet::default();
    let mut subspecies: HashSet<String> = HashSet::default();
    let mut strain: HashSet<String> = HashSet::default();

    let mut map = HashMap::default();
    let mut map_inv: HashMap<String, u32> = HashMap::default();
    let mut taxonomy_graph: HashMap<u32, u32> = HashMap::default(); //unidirectional graph connecting lower to higher order taxa

    //let mut species_lineage = HashMap::default();
    //let mut set: HashSet<String> = HashSet::default();
    //get sets for every taxonomic level
    for lineage in taxonomy_vector {
        let v: Vec<&str> = lineage.split(';').collect();
        //set.insert(v.join(";"));
        for i in (1..v.len()).rev() {
            match i {
                8 => strain.insert(v[..i + 1].join(";")),
                7 => subspecies.insert(v[..i + 1].join(";")),
                6 => species.insert(v[..i + 1].join(";")),
                5 => genus.insert(v[..i + 1].join(";")),
                4 => family.insert(v[..i + 1].join(";")),
                3 => order.insert(v[..i + 1].join(";")),
                2 => class.insert(v[..i + 1].join(";")),
                1 => phylum.insert(v[..i + 1].join(";")),
                _ => continue,
            };
        }
        kingdom.insert(v[0].to_string());
    }
    eprintln!("strain: {}, subspecies: {}, species: {}, genus: {}, family {}, order {}, class {}, phylum {}, kingdom {} ",
              strain.len(), subspecies.len(), species.len(), genus.len(), family.len(), order.len(), class.len(), phylum.len(), kingdom.len() );
    let mut counter: u32 = 1; //start with 1, reserve 0 for no hits/data
    for t in &kingdom {
        map_inv.insert(t.to_owned(), counter);
        map.insert(counter, t.to_owned());
        counter += 1;
    }
    for t in &phylum {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &class {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &order {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &family {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &genus {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &species {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &subspecies {
        map.insert(counter, t.to_owned());
        map_inv.insert(t.to_owned(), counter);
        counter += 1;
    }
    for t in &strain {
        map_inv.insert(t.to_owned(), counter);
        map.insert(counter, t.to_owned());
        counter += 1;
    }
    for lineage in taxonomy_vector {
        let v: Vec<&str> = lineage.split(';').collect();
        for i in (1..v.len()).rev() {
            match i {
                8 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                7 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                6 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                5 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                4 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                3 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                2 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                1 => taxonomy_graph
                    .insert(map_inv[&v[..i + 1].join(";")], map_inv[&v[..i].join(";")]),
                _ => continue,
            };
        }
    }
    let unclassified = (map.len() + 1) as u32; //maybe root?
    map.insert(0 as u32, "no hits".to_string());
    map_inv.insert("no hits".to_string(), 0 as u32);
    map.insert(unclassified, "unclassified".to_string());
    map_inv.insert("unclassified".to_string(), unclassified);
    (map, map_inv, taxonomy_graph)
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
