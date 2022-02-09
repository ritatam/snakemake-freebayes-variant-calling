import glob

configfile: "config.yaml"

SAMPLES = config["SRA"]     # Subset of samples to use as a list specified in config.yaml

def grep_input_to_merge(map_out):
    """Used in rule samtools_merge. Creates a string of absolute paths to polished bam files separated by space."""
    return ' '.join(glob.glob(f"{map_out}/*/*.dup_removed.bam"))



rule all:
    input:
        expand("{path}/{sample}/{ref}.{sample}.stat", path=config["map_out"], sample=SAMPLES, ref=config["bwa_ref_name"]),
        config["vcf_out"] + "/" + config["bwa_ref_name"] + ".merged.freebayes.vcf"




# Index the reference genome with bwa-mem2
rule bwamem2_index:
    output:
        index = multiext(config["bwa_ref_fa"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {config[env]} "
        " && echo $CONDA_PREFIX; "
        
        "bwa-mem2 index {config[bwa_ref_fa]}"


# Map paired-end reads to reference genome. Output bam sorted and indexed
rule bwamem2_map:
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:{sample}_lib1"
    input:
        index = multiext(config["bwa_ref_fa"], ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        r1 = config["trimmed_in"] + "/{sample}/{sample}_1.fastq.gz",
        r2 = config["trimmed_in"] + "/{sample}/{sample}_2.fastq.gz"
    output:
        bam = config["map_out"] + "/{sample}/" + config["bwa_ref_name"] + ".bwa_mem2.{sample}.bam",
    envmodules:
        "samtools/1.12"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {config[env]} "
        " && echo $CONDA_PREFIX; "
        
        "bwa-mem2 mem -R '{params.rg}' -t {config[bwa_mapT]} {config[bwa_ref_fa]} {input.r1} {input.r2} | samtools sort -O BAM -@ {config[bwa_sortT]} -o {output.bam} -"
        " && samtools index -@ {config[bwa_indexT]} {output.bam}"


# Collect statistics from bam
rule samtools_stat:
    input:
        bam = config["map_out"] + "/{sample}/" + config["bwa_ref_name"] + ".bwa_mem2.{sample}.bam"
    output:
        samstats = config["map_out"] + "/{sample}/" + config["bwa_ref_name"] + ".{sample}.stat"
    envmodules:
        "samtools/1.12"
    shell:
        "samtools stats -@ {config[sam_statT]} {input.bam} > {output.samstats}"


# Remove PCR duplicates from bam with Picard
rule picard_remove_duplicates:
    input:
        inbam = config["map_out"] + "/{sample}/" +  config["bwa_ref_name"] + ".bwa_mem2.{sample}.bam"
    output:
        outbam = config["map_out"] + "/{sample}/" + config["bwa_ref_name"] + ".bwa_mem2.{sample}.picard.dup_removed.bam",
        metrics = config["map_out"] + "/{sample}/" + config["bwa_ref_name"] + ".{sample}.picard.marked_dup_metrics.txt"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {config[env]} "
        " && echo $CONDA_PREFIX; "
        
        "picard MarkDuplicates -I {input.inbam} -O {output.outbam} -M {output.metrics} --REMOVE_DUPLICATES true"


# Merge polished bam from different samples into a single bam
rule samtools_merge:
    params:
        grep_input_to_merge(config["map_out"])
    input:
        expand(config["map_out"] + "/{sample}/" +  config["bwa_ref_name"] + ".bwa_mem2.{sample}.picard.dup_removed.bam", sample=SAMPLES)
    output:
        mergedbam = config["merged_out"] + "/" + config["bwa_ref_name"] + ".merged.bam"
    envmodules:
        "samtools/1.12"
    shell:
        "samtools merge -@ {config[sam_mergeT]} {output.mergedbam} {params}"


# Index the merged bam
rule samtools_merge_index:
    input:
        mergedbam = config["merged_out"] + "/" + config["bwa_ref_name"] + ".merged.bam"
    output:
        mergedbai = config["merged_out"] + "/" + config["bwa_ref_name"] + ".merged.bam.bai"
    envmodules:
        "samtools/1.12"
    shell:
        "samtools index -@ {config[sam_indexT]} {input.mergedbam}"


# Index the reference genome for freebayes
rule samtools_faidx:
    output:
        ref = config["bwa_ref_fa"] + ".fai"
    envmodules:
        "samtools/1.12"
    shell:
        "samtools faidx {config[bwa_ref_fa]}"


# Variant calling with freebayes-parallel
rule freebayes_parallel:
    input:
        fai = config["bwa_ref_fa"] + ".fai",
        mergedbai = config["merged_out"] + "/" + config["bwa_ref_name"] + ".merged.bam.bai",
        mergedbam = config["merged_out"] + "/" + config["bwa_ref_name"] + ".merged.bam"
    output:
        vcf = config["vcf_out"] + "/" + config["bwa_ref_name"] + ".merged.freebayes.vcf"
    envmodules:
        "parallel/20191022"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh "
        " && conda activate {config[env]} "
        " && echo $CONDA_PREFIX; "
        
        "freebayes-parallel <(fasta_generate_regions.py {input.fai} 100000) {config[freebayesT]} -f {config[bwa_ref_fa]} {input.mergedbam} > {output.vcf}"
