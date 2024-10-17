if config["call_snps"]:

    rule bcftools_call_snp:
        input:
            bam="seqnado_output/aligned/{sample}.bam",
            bai="seqnado_output/aligned/{sample}.bam.bai",
        output:
            vcf="seqnado_output/variant/bcftools/{sample}.vcf.gz",
        params:
            fasta=config["fasta"],
            faidx=config["fasta_index"],
        resources:
            mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
            runtime=lambda wildcards, attempt: f"{6 * 2 ** (attempt - 1)}h",
        threads: config["bcftools"]["threads"]
        log:
            "seqnado_output/logs/variant/bcftools/bcftools/{sample}.log",
        shell:
            "bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1"

    rule bcftools_split_multiallelic:
        input:
            vcf=rules.bcftools_call_snp.output.vcf,
        output:
            vcf="seqnado_output/variant/bcftools/{sample}.split.vcf.gz",
        resources:
            mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
            runtime=lambda wildcards, attempt: f"{6 * 2 ** (attempt - 1)}h",
        threads: config["bcftools"]["threads"]
        shell:"""
        bcftools norm --threads {threads} -Ou -m - {input.vcf} -Oz -o {output.vcf}
        tabix -f {output.vcf} > {output.vcf}.tbi
        """


    rule bcftools_filter_snp:
        input:
            vcf=rules.bcftools_split_multiallelic.output.vcf,
        output:
            vcf="seqnado_output/variant/bcftools/{sample}.filtered.vcf.gz",
        params:
            options=check_options(config["bcftools"]["filter"]),
        resources:
            mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
            runtime=lambda wildcards, attempt: f"{6 * 2 ** (attempt - 1)}h",
        threads: config["bcftools"]["threads"]
        shell: """
        bcftools view --threads {threads} {params.options} -Oz -o {output.vcf} {input.vcf} 
        tabix -f {output.vcf} > {output.vcf}.tbi
        """

    rule bcftools_annotate:
        input:
            vcf=rules.bcftools_call_snp.output.vcf,
            idx=rules.index_snp.output.vcf,
        output:
            vcf="seqnado_output/variant/bcftools/{sample}.anno.vcf.gz",
        params:
            dbsnp=config["snp_database"],
        threads: config["bcftools"]["threads"]
        shell:
            "bcftools annotate --threads {threads} -a {params.dbsnp} -c ID -Oz -o {output.vcf} {input.vcf}"

    rule bcftools_stats:
        input:
            vcf=rules.bcftools_annotate.output.vcf,
        output:
            stats="seqnado_output/variant/bcftools/{sample}.anno.stats.txt",
        params:
            fasta=config["fasta"],
        threads: config["bcftools"]["threads"]
        shell:
            "bcftools stats -F {params.fasta} -s - {input.vcf} > {output.stats}"

    rule bcftools_stats_plot:
        input:
            stats=rules.bcftools_stats.output.stats,
        output:
            summary="seqnado_output/variant/bcftools/{sample}/{sample}_summary.pdf",
        params:
            fasta=config["fasta"],
            out_dir="seqnado_output/variant/bcftools/{sample}",
        shell:
            "plot-vcfstats -p {params.out_dir} {input.stats}"