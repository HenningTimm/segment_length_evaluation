use clap::load_yaml;

use segment_length_analysis::hash_functions::{HlinParams, InvMultParams, Mmh3Params};
use segment_length_analysis::qgram_iterator::{Canonical, HashFunction};
use tab_hash;

fn main() {
    color_backtrace::install();

    let args_yml = load_yaml!("cl_args.yaml");
    let matches = clap::App::from_yaml(args_yml).get_matches();

    let fasta_path = matches.value_of("fasta-path").unwrap();
    let seg_len_path = matches.value_of("seg-len-path").unwrap();
    let minimizer_path = matches.value_of("minimizer-path").unwrap();
    let q = matches.value_of("q").unwrap().parse::<usize>().unwrap();
    let w = matches.value_of("w").unwrap().parse::<usize>().unwrap();
    let canonical;
    match matches.value_of("canonical").unwrap() {
        "non" => canonical = Canonical::No,
        "min" => canonical = Canonical::Min,
        "max" => canonical = Canonical::Max,
        _ => panic!(
            "Invalid value for canonical: {}",
            matches.value_of("canonical").unwrap()
        ),
    }
    let hash_function;
    match matches.value_of("mixing").unwrap() {
        "2bit" => hash_function = HashFunction::No,
        "fmix" => hash_function = HashFunction::Fmix,
        "swap" => hash_function = HashFunction::WordSwap,
        "hlin" => hash_function = HashFunction::Hlin(HlinParams::new()),
        "mmh3" => hash_function = HashFunction::Mmh3(Mmh3Params::new()),
        "tab64simple"  => hash_function = HashFunction::Tab64Simple(tab_hash::Tab64Simple::new()),
        "tab64twisted"  => hash_function = HashFunction::Tab64Twisted(tab_hash::Tab64Twisted::new()),
        "tab-reduced"  => hash_function = HashFunction::TabReduced(tab_hash::Tab64Twisted::new()),
        "inv-mult"  => hash_function = HashFunction::InvMult(InvMultParams::new()),
        _ => panic!("Invalid value for hash-mixing"),
    }
    let canonisize_qgrams = matches.is_present("canonisize-qgrams");
    let write_minimizers = matches.is_present("write-minimizers");
    segment_length_analysis::count_segment_lengths(
        fasta_path,
        seg_len_path,
        minimizer_path,
        q,
        w,
        canonical,
        hash_function,
        canonisize_qgrams,
        write_minimizers,
    );
}
