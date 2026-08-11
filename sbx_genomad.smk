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
            VIRUS_FP / "genomad" / "{sample}_genomad_virus_summary.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "genomad" / "{sample}_genomad_virus.fna",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "genomad" / "{sample}_genomad_virus_proteins.faa",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "genomad" / "{sample}_genomad_proteins.faa",
            sample=Samples.keys(),
        ),
        VIRUS_FP / "genomad" / "genomad_virus_summary.tsv",
        VIRUS_FP / "genomad" / "genomad_virus_genomes.fna",
        VIRUS_FP / "genomad" / "genomad_vOTU_repr_contig_map.tsv",
        VIRUS_FP / "genomad" / "genomad_vOTU_repr_contigs.fna",
        expand(
            VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.idxstats.tsv",
            sample=Samples.keys(),
        ),
        expand(
            VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.genomecov.tsv",
            sample=Samples.keys(),
        ),
        VIRUS_FP / "genomad" / "genomad_vOTU_mapped_read_counts.tsv",
        VIRUS_FP / "genomad" / "genomad_vOTU_taxonomy.tsv",


rule run_genomad:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
    output:
        summary=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_summary.tsv",
        vfna=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus.fna",
        vfaa=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_proteins.faa",
        faa=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_annotate" / "final.contigs_proteins.faa",
    benchmark:
        BENCHMARK_FP / "run_genomad_{sample}.tsv"
    log:
        LOG_FP / "run_genomad_{sample}.log",
    params:
        out_dir=str(VIRUS_FP / "genomad"),
        sample="{sample}",
        db_fp=Cfg["sbx_genomad"]["genomad_db"],
    conda:
        "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-genomad"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=1440,
        cpus_per_task=8,
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
        genomad end-to-end --cleanup {input.contigs} {params.sample} {params.db_fp} --threads {threads} >> {log} 2>&1
        """

rule genomad_output_format:
    input:
        summary=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_summary.tsv",
        vfna=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus.fna",
        vfaa=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_summary" / "final.contigs_virus_proteins.faa",
        faa=VIRUS_FP / "genomad" / "{sample}" / "final.contigs_annotate" / "final.contigs_proteins.faa",
    output:
        prefixed_summary=VIRUS_FP / "genomad" / "{sample}_genomad_virus_summary.tsv",
        prefixed_vfna=VIRUS_FP / "genomad" / "{sample}_genomad_virus.fna",
        prefixed_vfaa=VIRUS_FP / "genomad" / "{sample}_genomad_virus_proteins.faa",
        prefixed_faa=VIRUS_FP / "genomad" / "{sample}_genomad_proteins.faa",
    benchmark:
        BENCHMARK_FP / "genomad_fmt_{sample}.tsv"
    params:
        sample="{sample}",
    conda:
        "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-genomad"
    threads: 1
    resources:
        mem_mb=1000,
        runtime=60,
    shell:
        """
        # Add sample ID as new column in summary table and prefix sample ID to contig IDs
        awk -v new_contig_column="sample_contig" -v new_col_name="sample" -v new_col={params.sample} '
        BEGIN {{FS="\t"; OFS="\t"}} 
        NR==1 {{print new_contig_column, $0, new_col_name; next}}
        {{print new_col "_" $1, $0, new_col}}
        ' {input.summary} > {output.prefixed_summary} 

        # Prefix sample ID to virus fasta headers
        awk -v prefix={params.sample}_ '
        /^>/ {{sub(/^>/,">"prefix,$1)}}1
        ' {input.vfna} > {output.prefixed_vfna}

        # Prefix sample ID to virus proteins fasta headers
        awk -v prefix={params.sample}_ '
        /^>/ {{sub(/^>/,">"prefix,$1)}}1
        ' {input.vfaa} > {output.prefixed_vfaa} 
        
        # Prefix sample ID to all proteins fasta headers
        awk -v prefix={params.sample}_ '
        /^>/ {{sub(/^>/,">"prefix,$1)}}1
        ' {input.faa} > {output.prefixed_faa}
        """

rule genomad_output_merge:
    input:
        tables=expand(
            VIRUS_FP / "genomad" / "{sample}_genomad_virus_summary.tsv",
            sample=Samples.keys(),
        ),
        seqs=expand(
            VIRUS_FP / "genomad" / "{sample}_genomad_virus.fna",
            sample=Samples.keys(),
        ),
    output:
        tables=VIRUS_FP / "genomad" / "genomad_virus_summary.tsv",
        seqs=VIRUS_FP / "genomad" / "genomad_virus_genomes.fna",
    threads: 1
    conda:
       "envs/genomad_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-genomad"
    shell:
        """ 
        awk 'FNR==1 {{ if (NR==1) print; next }} 1' {input.tables} > {output.tables}
        cat {input.seqs} > {output.seqs}
        """

rule genomad_cluster:
    input:
        VIRUS_FP / "genomad" / "genomad_virus_genomes.fna"
    output:
        fltr=temp(VIRUS_FP / "genomad" / "genomad_vOTU_fltr.txt"),
        ani=VIRUS_FP / "genomad" / "genomad_virus_genomes_ani.tsv",
        aniids=VIRUS_FP / "genomad" / "genomad_virus_genomes_ani.ids.tsv",
        map=VIRUS_FP / "genomad" / "genomad_vOTU_repr_contig_map.tsv",
        size=VIRUS_FP / "genomad" / "genomad_vOTU_num_contigs.tsv",
        seqs=VIRUS_FP / "genomad" / "genomad_vOTU_repr_contigs.fna",
        bwaidx=VIRUS_FP / "genomad" / "genomad_vOTU_repr_contigs.fna.amb",
    threads: 16
    resources:
        mem_mb=96000,
        runtime=1440,
        cpus_per_task=16,
    log:
        LOG_FP / "genomad_cluster_vOTUs.log",
    conda:
        "envs/vclust_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_genomad:{SBX_GENOMAD_VERSION}-genomad"
    shell:
        """
        # Create a pre-alignment filter.
        vclust prefilter -i {input} -o {output.fltr} --min-ident 0.95 -t {threads} >> {log} 2>&1
        # Calculate ANI measures for genome pairs specified in the filter.
        vclust align -i {input} -o {output.ani} --filter {output.fltr} -t {threads} >> {log} 2>&1
        # Cluster contigs into vOTUs using the MIUVIG thresholds and the Leiden algorithm.
        # vclust cluster -i {output.ani} -o {output.map} --ids {output.aniids} --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 >> {log} 2>&1
        # Cluster for representative contigs
        vclust cluster -i {output.ani} -o {output.map} --out-repr --ids {output.aniids} --algorithm leiden --metric ani --ani 0.95 --qcov 0.85 >> {log} 2>&1
        # Get only the representative clusters and count the number of contigs in the cluster
        echo -e "cluster\tcontig_count"; awk -F'\t' 'NR>1 {{count[$2]++}} END {{for (c in count) print c "\t" count[c]}}' {output.map} | sort -k2,2nr > {output.size} 2>> {log}
        # Put all the sequences of representative clusters into one fasta file by extracting matching headers
        seqkit grep -f <(tail -n +2 {output.size} | cut -f1) {input} --threads {threads} > {output.seqs} 2>> {log}
        bwa index {output.seqs} >> {log} 2>&1
        """
    
rule genomad_read_mapping:
    input:
        seqs=VIRUS_FP / "genomad" / "genomad_vOTU_repr_contigs.fna",
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        bam=VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.bam",
        bai=VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.bam.bai",
        stats=VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.idxstats.tsv",
        cov=VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.genomecov.tsv",
    conda:
        "envs/vclust_env.yml"
    threads: 4
    resources:
        cpus_per_task=4,
    log:
        LOG_FP / "genomad_read_mapping_{sample}.log",
    shell:
        """
        bwa mem -t {threads} {input.seqs} {input.r1} {input.r2} | samtools view -u - | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam} {output.bai} >> {log} 2>&1 
        samtools idxstats {output.bam} | awk -v sampid={wildcards.sample} 'BEGIN{{print "ref\tseqlen\t" sampid "\tunmapped"}} {{print}}'> {output.stats} 2>> {log}
        bedtools genomecov -ibam {output.bam} -bga | awk 'BEGIN{{print "genome\tposition_start\tposition_end\tdepth"}} {{print}}'> {output.cov} 2>> {log}
        """

rule genomad_vOTU_table:
    input:
        stats=expand(
            VIRUS_FP / "genomad" / "{sample}_reads_to_vOTUs.idxstats.tsv",
            sample=Samples.keys(),
        ),
        summary=VIRUS_FP / "genomad" / "genomad_virus_summary.tsv",
    output:
        otutab=VIRUS_FP / "genomad" / "genomad_vOTU_mapped_read_counts.tsv",
        taxtab=VIRUS_FP / "genomad" / "genomad_vOTU_taxonomy.tsv",
    threads: 1
    log:
        LOG_FP / "genomad_vOTU_table.log",
    conda:
        "envs/genomad_env.yml"
    shell:
        """ 
        # grab mapped reads and reference length for vOTU table
        paste {input.stats} | awk '{{printf "%s\t%s", $1, $2; for(i=3;i<=NF;i+=4) printf "\t%s", $i; print ""}}' > {output.otutab} 2>> {log}
        # grab the vOTU representative sequence taxonomic assignments by contig ID and put them in their own file
        head -n1 {input.summary} > {output.taxtab}
        grep -F "$(awk '{{print $1 "\\t"}}' {output.otutab})" {input.summary} >> {output.taxtab} 2>> {log}
        """

