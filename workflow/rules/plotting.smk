########################################################################################################################
# Plot empirical vs. predicted segment length distributions
########################################################################################################################


rule plot_predicted_distributions:
    input:
        distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(50, 5000, 100)]
    output:
        single_path="../results/plots/distribution/prob_distr_single.pdf",
        heatmap_path="../results/plots/distribution/prob_distr_heatmap.pdf",
        lineplot_path="../results/plots/distribution/prob_distr_lineplot.pdf",
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/visualize_predicted_distributions.py"


rule plot_successive_differences_between_predictions:
    input:
        distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(50, 30000, 100)]
    output:
        errors_path="../results/plots/distribution/successive_errors.pdf",
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_successive_errors.py"


rule plot_errors_vs_empiric_distribution:
    input:
        small_distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(50, 5000, 1000)],
        large_distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(5000, 30001, 5000)],
        empiric = f"{tmp_prefix}/counts/synth/q=31_w=50_canon=non_hf=tab64twisted_genome=simulated_seg_lens.dat",
    output:
        errors_path="../results/plots/distribution/errors_vs_empiric.pdf",
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_empiric_errors.py"

# This is almost the same as above, but uses a special hash function tab-reduced
# which is a tabulation hash fucntion with a codomain size of 30011.
rule plot_errors_vs_empiric_distrubution_reduced:
    input:
        small_distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(50, 5000, 1000)],
        large_distribution_files=lambda wc: [f"../segdist-ng/distributions/prediction_w=50_k=150_C={C}.txt" for C in range(5000, 30001, 5000)],
        empiric = f"{tmp_prefix}/counts/synth/q=31_w=50_canon=non_hf=tab-reduced_genome=simulated_seg_lens.dat",
    output:
        errors_path="../results/plots/distribution/errors_vs_empiric_reduced.pdf",
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_empiric_errors.py"


rule plot_distribution_vs_reference_genomes:
    input:
        data_files=expand(f"{tmp_prefix}/counts/genomes/q={{{{q}}}}_w={{{{w}}}}_canon={{canonical}}_hf={{hf}}_genome={{{{genome}}}}_seg_lens.dat",
                           canonical=aps['canonical'],
                           hf=aps['hash_function'],
        ),
        predicted_distribution = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}.txt",
        predicted_error = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}_error.txt",
    output:
        landscape_pdf="../results/plots/genomes_landscape/q={q}_w={w}_genome={genome}_landscape_vs_prediction.pdf"
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_genome_landscape.py"


rule plot_prediction_vs_simulated_genomes:
    input:
        data_files=expand(f"{tmp_prefix}/counts/synth/q={{{{q}}}}_w={{{{w}}}}_canon={{canonical}}_hf={{hf}}_genome=simulated_seg_lens.dat",
                           canonical=aps['canonical'],
                           hf=aps['hash_function'],
        ),
        predicted_distribution = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}.txt",
        predicted_error = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}_error.txt",
    output:
        landscape_pdf="../results/plots/simulated_landscape/q={q}_w={w}_genome=simulated_landscape_vs_prediction.pdf"
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_genome_landscape.py"


rule plot_reduced_vs_simulated:
    input:
        data_files=expand(f"{tmp_prefix}/counts/synth/q={{{{q}}}}_w={{{{w}}}}_canon={{canonical}}_hf={{hf}}_genome=simulated_seg_lens.dat",
                           canonical=aps['canonical'],
                           hf=aps['hash_function_reduced'],
        ),
        predicted_distribution = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}.txt",
        predicted_error = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}_error.txt",
    output:
        landscape_pdf="../results/plots/simulated_landscape_reduced/q={q}_w={w}_genome=simulated_reduced_landscape_vs_prediction.pdf"
    conda:
        "../envs/plot_env.yaml"
    params:
        fontscale=2
    script:
        "../scripts/plot_genome_landscape.py"



rule plot_reference_genomes_comparison:
    input:
        data_files=expand(f"{tmp_prefix}/counts/genomes/q={{{{q}}}}_w={{{{w}}}}_canon={{{{canonical}}}}_hf={{hf}}_genome={{genome}}_seg_lens.dat",
                          genome=aps["genome_comp"],
                          hf=aps['genome_comp_hfs'],
        ),
        predicted_distribution = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}.txt",
        predicted_error = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}_error.txt",
    output:
        genome_comp_pdf="../results/plots/genomes_comparison/q={q}_w={w}_canon={canonical}_genome_comparison.pdf"
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_genome_comparison.py"


########################################################################################################################
# Plot different GC values
########################################################################################################################

rule plot_gc_overview:
    input:
        data_files=expand(f"{tmp_prefix}/counts/gc_content_analysis/q={{{{q}}}}_w={{{{w}}}}_canon={{canonical}}_hf={{hf}}_gc={{gc}}_seg_lens.dat",
               gc=config['synthetic_genomes']['gc_content_range'],
               canonical=aps['canonical'],
               hf=aps['hash_function'],
        ),
        predicted_distribution = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}.txt",
        predicted_error = lambda wc: f"../segdist-ng/distributions/prediction_w={wc.w}_k={int(wc.w)*3}_C={aps['ref_C']}_error.txt",        
    output:
        gc_landscape_pdf="../results/plots/gc_content_analysis/overview_w={w}_q={q}.pdf"
    params:
        combinations=[
            ("2bit", "non"), ("2bit", "min"), ("2bit", "max"),
            ("tab64twisted", "non"), ("tab64twisted", "min"), ("tab64twisted", "max"),
            ("mmh3", "min"),
        ]
    conda:
        "../envs/plot_env.yaml"
    script:
        "../scripts/plot_gc_facets.py"
