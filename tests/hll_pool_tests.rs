use sepia::hll_pool::{get_hll, put_hll, init_pool};

#[test]
fn test_hll_pool_basic() {
    // ensure pool has some capacity
    init_pool(4);
    let mut h = get_hll();
    // initially empty
    assert_eq!(h.len().trunc() as usize, 0);
    h.insert(&42u64);
    // now non-empty
    assert!(h.len().trunc() as usize >= 1);
    put_hll(h);
    // taking again should give an empty HLL (pool returns a reset instance)
    let h2 = get_hll();
    assert_eq!(h2.len().trunc() as usize, 0);
}
