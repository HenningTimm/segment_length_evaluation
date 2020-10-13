# Genome Simulation
This Rust program simulates a single chromosome with `l` bases and a GC-content of `g` in FASTA format.
The simulated sequence is written to stdout.
Usage:

```
cargo run --release -- -l 100 -g 0.75 > 100bases_gc075.fasta
```

The above command generates one sequence with 100 bases, 75% of which are `G` or `C`.

