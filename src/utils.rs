use flate2::{bufread, write, Compression};
use std::fs::File;
use std::io;
use std::path::Path;

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


pub fn get_tsv_reader<P: AsRef<Path>>(path: P) -> Result<csv::Reader<Box<dyn io::Read>>, io::Error> {
    let fh = get_reader_handle(path.as_ref().to_str().expect("Not UTF-8 file name"))?;
    let reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(fh);
    Ok(reader)
}

pub fn get_tsv_headers<R: io::Read>(
    reader: &mut csv::Reader<R>,
) -> Result<csv::ByteRecord, csv::Error> {
    let headers = reader.headers().unwrap().as_byte_record().to_owned();
    Ok(headers)
}