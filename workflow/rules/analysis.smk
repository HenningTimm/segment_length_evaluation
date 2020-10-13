####################################################################################################
# Compute segment length counts
####################################################################################################

rule analyze_segment_lengths_for_simulated_genome:
    input:
        "seg_len_analyzer",
        fasta_file = f"{simulation_prefix}/sg_{{nr}}.fasta",
    output:
        seg_lens=f"{tmp_prefix}/counts/synth/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_nr={{nr}}_seg_lens.dat",
        minimizers=f"{tmp_prefix}/counts/synth/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_nr={{nr}}_minimizers.dat",
    conda:
        "../envs/rust_env.yaml"
    shell:
        "./seg_len_analyzer {input.fasta_file} {output.seg_lens} {output.minimizers} "
        "-q {wildcards.q} -w {wildcards.w} --canonical {wildcards.canonical} "
        "--hash-mixing {wildcards.hf} --canonisize-qgrams"


rule analyze_segment_lengths_for_gc_content:
    input:
        "seg_len_analyzer",
        fasta_file = f"{simulation_prefix}/gc_content_analysis/synth_gc={{gc}}.fasta",
    output:
        seg_lens=f"{tmp_prefix}/counts/gc_content_analysis/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_gc={{gc}}_seg_lens.dat",
        minimizers=f"{tmp_prefix}/counts/gc_content_analysis/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_gc={{gc}}_minimizers.dat",
    conda:
        "../envs/rust_env.yaml"
    shell:
        "./seg_len_analyzer {input.fasta_file} {output.seg_lens} {output.minimizers} "
        "-q {wildcards.q} -w {wildcards.w} --canonical {wildcards.canonical} "
        "--hash-mixing {wildcards.hf} --canonisize-qgrams"


rule analyze_segment_lengths_for_reference_genome:
    input:
        "seg_len_analyzer",
        fasta_file=lambda w: g[w.genome]["path"],
    output:
        seg_lens=f"{tmp_prefix}/counts/genomes/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_genome={{genome}}_seg_lens.dat",
        minimizers=f"{tmp_prefix}/counts/genomes/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_genome={{genome}}_minimizers.dat",
    wildcard_constraints:
        q="\d+",
        w="\d+",
    conda:
        "../envs/rust_env.yaml"
    shell:
        "./seg_len_analyzer {input.fasta_file} {output.seg_lens} {output.minimizers} "
        "-q {wildcards.q} -w {wildcards.w} --canonical {wildcards.canonical} "
        "--hash-mixing {wildcards.hf} --canonisize-qgrams"



####################################################################################################
# Aggregate counts for plotting
####################################################################################################


rule aggregate_counts_for_simulated_genomes:
    input:
        expand(f"{tmp_prefix}/counts/synth/q={{{{q}}}}_w={{{{w}}}}_canon={{{{canonical}}}}_hf={{{{hf}}}}_nr={{nr}}_seg_lens.dat",
               nr=range(0, config['synthetic_genomes']['nr']),
        )
    output:
        f"{tmp_prefix}/counts/synth/q={{q}}_w={{w}}_canon={{canonical}}_hf={{hf}}_genome=simulated_seg_lens.dat"
    conda:
        "../envs/plot_env.yaml"
    shell:
        "python scripts/aggregate_counts.py {input} --output_path {output}"
