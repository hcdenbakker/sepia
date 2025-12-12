use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;

/// Start a writer thread that receives `Vec<u8>` chunks and writes them
/// to `filename`. If `compress` is true the output is gzipped.
pub fn start_writer_thread(
    filename: String,
    compress: bool,
) -> (Sender<Vec<u8>>, thread::JoinHandle<()>) {
    let (tx, rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = mpsc::channel();
    let handle = thread::spawn(move || {
        let file = File::create(filename).expect("could not create outfile!");
        if compress {
            let mut enc = GzEncoder::new(BufWriter::new(file), Compression::default());
            while let Ok(chunk) = rx.recv() {
                let _ = enc.write_all(&chunk);
            }
            let _ = enc.finish();
        } else {
            let mut w = BufWriter::new(file);
            while let Ok(chunk) = rx.recv() {
                let _ = w.write_all(&chunk);
            }
            let _ = w.flush();
        }
    });
    (tx, handle)
}
