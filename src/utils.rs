use flate2::{bufread, write, Compression};
use std::fs::File;
use std::io;

pub fn get_reader_handle(path: &str) -> Result<Box<dyn io::Read>, io::Error> {
    if path.ends_with(".gz") {
        let f = File::open(path).unwrap();
        Ok(Box::new(bufread::GzDecoder::new(io::BufReader::new(f))))
    } else {
        Ok(Box::new(File::open(path).unwrap()))
    }
}

pub fn get_writer_handle(path: &str) -> Result<Box<dyn io::Write>, io::Error> {
    if path.ends_with(".gz") {
        let f = File::create(path).unwrap();
        Ok(Box::new(write::GzEncoder::new(io::BufWriter::new(f), Compression::default())))
    } else {
        Ok(Box::new(File::create(path).unwrap()))
    }
}
