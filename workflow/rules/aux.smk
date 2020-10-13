##########################################################################################
# Build scripts for Rust code
##########################################################################################

# Build the executable used to analyze segment lengths
rule build_segment_length_analyzer:
    input:
        source_files = glob.glob("../segment_length_analysis/src/*.rs"),
    output:
        executable = "../segment_length_analysis/target/release/segment_length_analysis",
        symlink = "seg_len_analyzer",
    conda:
        "../envs/rust_env.yaml"
    shell:
        """
        echo {input.source_files}
        cd ../segment_length_analysis
        cargo build --release
        cd -
        ln -s {output.executable} {output.symlink}
        """

# Build the executable used to simulate genomes
rule build_genome_simulator:
    input:
        source_files = glob.glob("../genome_simulation/src/*.rs"),
    output:
        executable = "../genome_simulation/target/release/genome_simulation",
        symlink = "simulate_genome",
    conda:
        "../envs/rust_env.yaml"
    shell:
        """
        echo {input.source_files}
        cd ../genome_simulation
        cargo build --release
        cd -
        ln -s {output.executable} {output.symlink}
        """


##########################################################################################
# Predicted distributions
##########################################################################################

# Generate a predicted segment length distribution for a given hash function codomain size C
# winnowing window size w, and maximum segment length k.
rule generate_distribution_file:
    output:
        distr_file = "../segdist-ng/distributions/prediction_w={w}_k={k}_C={C}.txt",
        error_file = "../segdist-ng/distributions/prediction_w={w}_k={k}_C={C}_error.txt",
    conda:
        "../envs/plot_env.yaml"
    shell:
        """
        python ../segdist/segdist.py -w {wildcards.w} -K {wildcards.k} -C {wildcards.C} > {output.distr_file} 2> {output.error_file}
        """
