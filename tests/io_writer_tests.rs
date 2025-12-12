use sepia::io_writer::start_writer_thread;
use std::fs;
use flate2::read::GzDecoder;
use std::io::Read;

#[test]
fn test_io_writer_uncompressed() {
    let fname = "/tmp/sepia_test_io.txt".to_string();
    let (tx, h) = start_writer_thread(fname.clone(), false);
    let _ = tx.send(b"hello\n".to_vec());
    let _ = tx.send(b"world\n".to_vec());
    drop(tx);
    let _ = h.join();
    let content = fs::read_to_string(fname).expect("read failed");
    assert!(content.contains("hello"));
    assert!(content.contains("world"));
}

#[test]
fn test_io_writer_compressed() {
    let fname = "/tmp/sepia_test_io.gz".to_string();
    let (tx, h) = start_writer_thread(fname.clone(), true);
    let _ = tx.send(b"line1\n".to_vec());
    let _ = tx.send(b"line2\n".to_vec());
    drop(tx);
    let _ = h.join();
    let data = fs::read(fname).expect("read gz failed");
    let mut d = GzDecoder::new(&data[..]);
    let mut s = String::new();
    d.read_to_string(&mut s).expect("decompress failed");
    assert!(s.contains("line1"));
    assert!(s.contains("line2"));
}
