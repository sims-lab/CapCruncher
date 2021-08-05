pub mod count_fragments;



fn main(){

    let fragments = "/stopgap/milne_group/asmith/tiledc/capcruncher_analysis/reporters/SRR10095377.alpha-globin.tsv.gz";

    println!("Started loading tsv");
    let df = count_fragments::load_lazy_dataframe(fragments).expect("Error reading tsv");
    println!("Loaded tsv");
    println!("Started counting");
    let counts = count_fragments::count_fragments(df).expect("Error counting");
    println!("Writing to file");
    count_fragments::to_file("RUST_COUNTS.tsv", counts).expect("Error writing to file");


}