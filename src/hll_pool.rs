use std::cell::RefCell;
use hyperloglog::HyperLogLog;

thread_local! {
    static HLL_POOL: RefCell<Vec<HyperLogLog<u64>>> = RefCell::new(Vec::new());
}

/// Prefill the per-thread pool with `size` empty HLLs.
pub fn init_pool(size: usize) {
    HLL_POOL.with(|p| {
        let mut pool = p.borrow_mut();
        while pool.len() < size {
            pool.push(HyperLogLog::new(0.001));
        }
    });
}

/// Pop an HLL from the pool or create a new empty one.
pub fn get_hll() -> HyperLogLog<u64> {
    HLL_POOL.with(|p| p.borrow_mut().pop().unwrap_or_else(|| HyperLogLog::new(0.001)))
}

/// Return an HLL to the pool. The HLL's contents are cleared (replaced)
/// to avoid leaking previous data into future users.
pub fn put_hll(_h: HyperLogLog<u64>) {
    // Reset by creating a fresh, empty HLL and pushing to pool.
    HLL_POOL.with(|p| p.borrow_mut().push(HyperLogLog::new(0.001)));
}
