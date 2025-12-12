//https://stackoverflow.com/questions/29037033/how-to-slice-a-large-veci32-as-u8
extern crate byteorder;

use byteorder::{LittleEndian, WriteBytesExt};

use std::fs::File;
use std::io::{Read, Write};
use std::{mem, slice};

use std::io::BufWriter;
use memmap2::Mmap;

fn as_u8_slice(v: &[u32]) -> &[u8] {
    let element_size = mem::size_of::<u32>();
    unsafe { slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

fn from_u8(v: Vec<u8>) -> Vec<u16> {
    let data = v.as_ptr();
    let len = v.len();
    let capacity = v.capacity();
    let element_size = mem::size_of::<u16>();

    // Make sure we have a proper amount of capacity (may be overkill)
    assert_eq!(capacity % element_size, 0);
    // Make sure we are going to read a full chunk of stuff
    assert_eq!(len % element_size, 0);

    unsafe {
        // Don't allow the current vector to be dropped
        // (which would invalidate the memory)
        mem::forget(v);

        Vec::from_raw_parts(
            data as *mut u16,
            len / element_size,
            capacity / element_size,
        )
    }
}

fn from_u8_to_u32(v: Vec<u8>) -> Vec<u32> {
    let data = v.as_ptr();
    let len = v.len();
    let capacity = v.capacity();
    let element_size = mem::size_of::<u32>();

    // Make sure we have a proper amount of capacity (may be overkill)
    assert_eq!(capacity % element_size, 0);
    // Make sure we are going to read a full chunk of stuff
    assert_eq!(len % element_size, 0);

    unsafe {
        // Don't allow the current vector to be dropped
        // (which would invalidate the memory)
        mem::forget(v);

        Vec::from_raw_parts(
            data as *mut u32,
            len / element_size,
            capacity / element_size,
        )
    }
}

pub fn do_write_u32(filename: &str, v: &[u32]) {
    let mut f = File::create(filename).unwrap();
    f.write_all(as_u8_slice(v)).unwrap();
}

pub fn do_read(filename: &str) -> Vec<u16> {
    let mut f = File::open(filename).unwrap();
    let mut bytes = Vec::new();

    f.read_to_end(&mut bytes).unwrap();

    from_u8(bytes)
}

pub fn do_read_u32(filename: &str) -> Vec<u32> {
    let mut f = File::open(filename).unwrap();
    let mut bytes = Vec::new();

    f.read_to_end(&mut bytes).unwrap();

    from_u8_to_u32(bytes)
}

pub fn do_read_u32_mmap(filename: &str) -> Mmap {
    let file = File::open(filename).unwrap();
    unsafe { Mmap::map(&file).unwrap() }
}

pub fn mmap_as_u32_slice(mmap: &Mmap) -> &[u32] {
    let data = mmap.as_ptr();
    let len = mmap.len();
    let element_size = mem::size_of::<u32>();
    assert_eq!(len % element_size, 0);
    unsafe { slice::from_raw_parts(data as *const u32, len / element_size) }
}

pub fn do_write_little_endian(filename: &str, v: &[u16]) {
    let mut f = File::create(filename).expect("problems creating db file");
    let mut result: Vec<u8> = Vec::new();
    for &n in v {
        let _ = result.write_u16::<LittleEndian>(n);
    }
    f.write_all(&result[..]).unwrap();
}

pub fn do_write_little_endian_u32(filename: &str, v: &[u32]) {
    let mut f = BufWriter::new(File::create(filename).expect("problems creating db file"));
    let mut result: Vec<u8> = Vec::new();
    for &n in v {
        let _ = result.write_u32::<LittleEndian>(n);
    }
    f.write_all(&result[..]).unwrap();
}
/*
let slice_u16: &[u16] = &*vec![1, 2, 3, 4, 5, 6];
    println!("u16s: {:?}", slice_u16);
    for &n in v {
        let _ = result.write_u16::<LittleEndian>(n);
    }

    let mut result: Vec<u8> = Vec::new();
    for &n in slice_u16 {
        let _ = result.write_u16::<LittleEndian>(n);
    }



fn main() {
    let v = vec![42; 10];
    do_write("vector.dump", &v);
    let v2 = do_read("vector.dump");

    assert_eq!(v, v2);
    println!("{:?}", v2)
}*/
