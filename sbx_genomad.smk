try:
    SBX_GENOMAD_VERSION = get_ext_version("sbx_genomad")
except (NameError, ValueError):
    # For backwards compatibility with older versions of Sunbeam
    SBX_GENOMAD_VERSION = "0.0.0"
VIRUS_FP = output_subdir(Cfg, "virus")


def get_extension_path() -> Path:
    return Path(__file__).parent.parent.resolve()
def get_rules_path() -> Path:
    return Path(__file__).resolve()

rule all_genomad:
    input:
        expand(
            VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_summary.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "genomad" / "{sample}_final.contigs_summary.tsv",
            sample=Samples.keys(),
        ),
        VIRUS_FP / "genomad" / "genomad_all_contigs_virus_summary.tsv",

rule genomad:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        summary=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_summary.tsv",
    benchmark:
        BENCHMARK_FP / "genomad_{sample}.tsv"
    log:
        LOG_FP / "genomad_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "genomad"),
        sample="{sample}",
        db_fp=Cfg["sbx_genomad"]["genomad_db"],
    conda:
        "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-cenote-taker"
    resources:
        mem_mb=96000,
        runtime=1440,
        threads=8
    shell:
        """
        SAMPLE={params.sample}
        if [ -s {input.contigs} ]; then
            echo "Contigs file exists and is not empty" >> {log}
        else
            echo "Contigs file is empty" >> {log}
            touch {output.summary}
            exit 0
        fi

        if [ ! -d {params.db_fp} ] || [ ! "$(ls -A {params.db_fp})" ]; then
            echo "geNomad database path {params.db_fp} is missing or empty" >> {log}
            exit 1
        fi

        cd {params.out_dir}
        genomad end-to-end --cleanup {input.contigs} {params.sample} {params.db_fp} --threads {resources.threads} >> {log} 2>&1
        """

rule genomad_fmt:
    input:
        summary=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_summary.tsv",
    output:
        prefixed_summary=VIRUS_FP / "genomad" / "{sample}_final.contigs_summary.tsv",
    benchmark:
        BENCHMARK_FP / "genomad_fmt_{sample}.tsv"
    params:
        sample="{sample}",
    conda:
        "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-cenote-taker"
    resources:
        mem_mb=1000,
        runtime=60,
        threads=1
    shell:
        """
        awk -v new_col_name="sample" -v new_col={params.sample} '
        BEGIN {{FS="\t"; OFS="\t"}} 
        NR==1 {{print $0, new_col_name; next}}
        {{print $0, new_col}}
        ' {input.summary} > {output.prefixed_summary} 
        """

rule genomad_summarize:
    input:
        expand(
            VIRUS_FP / "genomad" / "{sample}_final.contigs_summary.tsv",
            sample=Samples.keys(),
        ),
    output:
        VIRUS_FP / "genomad" / "genomad_all_contigs_virus_summary.tsv",
    conda:
       "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-cenote-taker"
    shell:
        """ 
        awk 'FNR==1 {{ if (NR==1) print; next }} 1' {input} > {output}
        """
