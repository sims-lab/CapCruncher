use itertools::Itertools;
use polars::prelude::*;
use std::collections::{BTreeMap};
use std::path::{Path};
use std::fs::File;
use std::io;
use std::io::Write;

pub fn load_dataframe(file: &str) -> Result<DataFrame> {
    let df = CsvReader::from_path(file)?
        .infer_schema(Some(10000))
        .with_delimiter(b'\t')
        .has_header(true)
        .finish()?;
    Ok(df)
}

pub fn load_lazy_dataframe(file: &str) -> Result<LazyFrame> {
    let df = LazyCsvReader::new(file.to_string())
        .with_delimiter(b'\t')
        .has_header(true)
        .finish();
    Ok(df)
}

pub fn count_fragments(df: LazyFrame) -> Result<BTreeMap<(i64, i64), usize>> {
    
    
    let groupby = df.groupby("parent_read");
    let groups = groupby.get_groups();
    let mut fragment_counts = BTreeMap::new();
    let restriction_fragments: Vec<_> = df
        .column("restriction_fragment")?
        .i64()
        .unwrap()
        .into_no_null_iter()
        .collect();

    for (_first_index, grouped_indicies) in groups {
        let fragments_for_combining = grouped_indicies
            .iter()
            .map(|ind| restriction_fragments[*ind as usize]);

        for mut fragments in fragments_for_combining.into_iter().combinations(2) {
            
            fragments.sort_unstable();
            let fragment1 = fragments[0];
            let fragment2 = fragments[1];

            let count = fragment_counts.entry((fragment1, fragment2)).or_insert(0);
            *count += 1;
        }
    }

    Ok(fragment_counts)
}

pub fn to_file<P: AsRef<Path>>(path: P, counts: BTreeMap<(i64, i64), usize>) -> std::result::Result<(), io::Error> {

    let outfile = File::create(path)?;
    let mut writer = io::BufWriter::new(outfile);

    writer.write_all(b"bin1_id\tbin2_id\tcount\n")?;
    for ((b1, b2), count) in counts.iter(){

        writer.write_all(format!("{}\t{}\t{}\n", b1, b2, count).as_bytes())?;
    

    } 

    Ok(())





}
