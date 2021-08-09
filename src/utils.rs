use flate2::{bufread, write, Compression};
use itertools::Itertools;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;

pub fn get_reader_handle(path: &str) -> Box<dyn io::Read> {
    if path.ends_with(".gz") {
        let f = File::open(path).unwrap();
        Box::new(bufread::GzDecoder::new(io::BufReader::new(f)))
    } else {
        Box::new(File::open(path).unwrap())
    }
}

pub fn get_writer_handle(path: &str) -> Box<dyn io::Write> {
    if path.ends_with(".gz") {
        let f = File::open(path).unwrap();
        Box::new(write::GzEncoder::new(io::BufWriter::new(f), Compression::default()))
    } else {
        Box::new(File::create(path).unwrap())
    }
}
