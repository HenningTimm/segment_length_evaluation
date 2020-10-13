import glob
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

##########################################################################################
# Genome simulation
##########################################################################################

rule simulate_default_genome:
    input:
        "simulate_genome"
    output:
        f"{simulation_prefix}/sg_{{nr}}.fasta"
    params:
        length = config['synthetic_genomes']['length'],
        gc_content = config['synthetic_genomes']['default_gc_content'],
    conda:
        "../envs/rust_env.yaml"
    shell:
        "./simulate_genome -l {params.length} --gc-content {params.gc_content} > {output[0]}"


rule simulate_gc_genome:
    input:
        "simulate_genome"
    output:
        f"{simulation_prefix}/gc_content_analysis/synth_gc={{gc}}.fasta"
    params:
        length = config['synthetic_genomes']['length'],
    conda:
        "../envs/rust_env.yaml"
    shell:
        "./simulate_genome -l {params.length} --gc-content {wildcards.gc} > {output[0]}"


##########################################################################################
# Genome Download
##########################################################################################

# NOTE: Download paths are specified in the file genomes.tsv

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()

rule get_mxanthus_genome:
    input:
        FTP.remote(
            g["mxanthus"]["download_link"],
            )
    output:
        g["mxanthus"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"


rule get_pfalciparum_genome:
    input:
        FTP.remote(
            g["pfalciparum"]["download_link"],
            )
    output:
        g["pfalciparum"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

rule get_hops_genome:
    input:
        HTTP.remote(
            g["hops"]["download_link"],
            insecure=True,  # Hopbase download is HTTP, not HTTPS
        )
    output:
        g["hops"]["path"],
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

rule get_human_genome:
    input:
        FTP.remote(
            g["hg38"]["download_link"],
        )
    output:
        g["hg38"]["path"],
    shell:
        "mkdir -p genomes && zcat {input} > {output}"
