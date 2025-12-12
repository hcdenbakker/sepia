use hashbrown::HashMap;

use sepia::taxonomy_u32;
use sepia::search_bits_sepia::poll_lineages_vec_u32;

// Build a simple taxonomy graph:
// root(1) -> A(2) -> A1(3)
//                \-> A2(4)
//         -> B(5) -> B1(6)

fn make_graph() -> (HashMap<u32,u32>, HashMap<u32,String>) {
    let mut graph = HashMap::new();
    // parent links: child -> parent
    graph.insert(2, 1);
    graph.insert(3, 2);
    graph.insert(4, 2);
    graph.insert(5, 1);
    graph.insert(6, 5);

    // taxonomy strings (semi-colon separated lineage strings)
    let mut taxonomy = HashMap::new();
    taxonomy.insert(1, "root".to_string());
    taxonomy.insert(2, "root;A".to_string());
    taxonomy.insert(3, "root;A;A1".to_string());
    taxonomy.insert(4, "root;A;A2".to_string());
    taxonomy.insert(5, "root;B".to_string());
    taxonomy.insert(6, "root;B;B1".to_string());
    taxonomy.insert(0, "no hits".to_string());

    (graph, taxonomy)
}

#[test]
fn test_find_lca_u32_basic() {
    let (graph, taxonomy) = make_graph();
    // unclassified token as used in code
    let unclassified = (taxonomy.len() + 1) as u32;

    // LCA of 3 (A1) and 4 (A2) should be 2 (A)
    let lca_a = taxonomy_u32::find_lca_u32(3, 4, &graph, unclassified);
    assert_eq!(lca_a, 2);

    // LCA of 3 (A1) and 6 (B1) should be root (1)
    let lca_root = taxonomy_u32::find_lca_u32(3, 6, &graph, unclassified);
    assert_eq!(lca_root, 1);

    // LCA of parent-child (2 and 3) -> 2
    let lca_parent = taxonomy_u32::find_lca_u32(2, 3, &graph, unclassified);
    assert_eq!(lca_parent, 2);
}

#[test]
fn test_poll_lineages_vec_u32_behavior() {
    let (graph, taxonomy) = make_graph();
    // Create a small report: (taxon, count)
    // Suppose observed counts: A1:5, A2:3, B1:2
    let report = vec![(3u32, 5usize), (4u32, 3usize), (6u32, 2usize)];

    // Call poll_lineages_vec_u32 and check winner and score
    let (winner, score, k) = poll_lineages_vec_u32(&report, 31, &graph, &taxonomy);

    // The aggregated scores per candidate (by string containment semantics):
    // For 3 (A1): contains A1,A => adds counts of A1 and A2 (5+3)=8 (doesn't contain B1)
    // For 4 (A2): contains A2,A => adds counts of A1 and A2 = 8
    // For 6 (B1): contains B1,B => adds counts of B1 =2
    // So top score ties between 3 and 4 -> LCA(3,4)=2
    assert_eq!(winner, 2);
    assert_eq!(score, 8);
    assert_eq!(k, 31);
}
