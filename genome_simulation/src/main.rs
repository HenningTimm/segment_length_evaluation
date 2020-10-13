#[macro_use]
extern crate clap;
extern crate rand;

use rand::distributions::{Weighted, WeightedChoice, Distribution};

fn main() {
    let args_yml = load_yaml!("cl_args.yaml");
    let matches = clap::App::from_yaml(args_yml).get_matches();
    let length = matches
        .value_of("length")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let gc_content = matches
        .value_of("gc_content")
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let gc_weight;
    let at_weight;
    if gc_content == 0.0 {
        gc_weight = 0.0;
        at_weight = 1.0;
    } else if gc_content == 1.0 {
        gc_weight = 1.0;
        at_weight = 0.0;        
    } else {
        // Map both weights to positive integers (required by weighted choice)
        let min = gc_content.min(1.0 - gc_content);
        gc_weight = gc_content / min;
        at_weight = (1.0 - gc_content) / min;
    }
    // Make sure the above calculation hold for untested combinations of gc_content
    // eprintln!("gc_weight: {}, at_weight: {}", gc_weight, at_weight);
    assert_eq!(gc_weight / at_weight, gc_content / (1.0 - gc_content));

    // Initialize RNG and weights to select bases
    let mut weighted_nts = vec![
        Weighted { weight: at_weight as u32, item: 65_u8 }, // A
        Weighted { weight: gc_weight as u32, item: 67_u8 }, // C
        Weighted { weight: gc_weight as u32, item: 71_u8 }, // G
        Weighted { weight: at_weight as u32, item: 84_u8 }, // T
    ];
    let nt_choice = WeightedChoice::new(&mut weighted_nts);
    let mut rng = rand::rngs::OsRng::new().unwrap();
    
    // Simulate bases
    let mut seq = Vec::with_capacity(length);
    for i in 0..length {
        if i % 10000 == 0 {
            eprint!("\r{:>12} / {:>12}", i, length)
        }
        seq.push(nt_choice.sample(&mut rng));
    }

    // Assemble header
    let name = format!(">Random_GC={}_len={}", gc_content, length);

    // Write output
    println!("{}", name);
    println!("{}", std::str::from_utf8(&seq).unwrap());
}
