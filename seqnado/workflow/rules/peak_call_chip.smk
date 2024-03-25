from typing import Literal
from seqnado.helpers import check_options
import re


def get_lanceotron_threshold(wildcards):
    options = config["lanceotron"]["callpeak"]
    threshold_pattern = re.compile(r"\-c\s+(\d+.?\d*)")
    threshold = threshold_pattern.search(options).group(1)
    return threshold


def get_control_bam(wildcards):
    exp = DESIGN.query(sample_name=wildcards.sample, ip=wildcards.treatment)
    if exp.control:
        control = f"seqnado_output/aligned/{wildcards.sample}_{exp.control}.bam"
    else:
        control = ""
    return control


def get_control_tag(wildcards):
    exp = DESIGN.query(sample_name=wildcards.sample, ip=wildcards.treatment)
    if not exp.control:
        control = ""
    else:
        control = f"seqnado_output/tag_dirs/{wildcards.sample}_{exp.control}"
    return control


def get_control_bigwig(wildcards):
    exp = DESIGN.query(sample_name=wildcards.sample, ip=wildcards.treatment)
    if not exp.control:
        control = ""
    else:
        control = f"seqnado_output/bigwigs/deeptools/unscaled/{wildcards.sample}_{exp.control}.bigWig"
    return control


rule macs2_with_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=get_control_bam,
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem="2GB",
        runtime="2h",
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """


rule macs2_no_input:
    input:
        treatment="seqnado_output/aligned/{sample}_{treatment}.bam",
        control=lambda wc: "UNDEFINED" if get_control_bam else [], 
    output:
        peaks="seqnado_output/peaks/macs/{sample}_{treatment}.bed",
    params:
        options=check_options(config["macs"]["callpeak"]),
        narrow=lambda wc, output: output.peaks.replace(".bed", "_peaks.narrowPeak"),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    resources:
        mem="2GB",
        runtime="2h",
    log:
        "seqnado_output/logs/macs/{sample}_{treatment}.log",
    shell:
        """
        macs2 callpeak -t {input.treatment} -n {params.basename} -f BAMPE {params.options} > {log} 2>&1 &&
        cat {params.narrow} | cut -f 1-3 > {output.peaks}
        """


rule homer_with_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=get_control_tag,
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem="4GB",
        runtime="2h",
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp  -i {input.control} > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule homer_no_input:
    input:
        treatment="seqnado_output/tag_dirs/{sample}_{treatment}",
        control=lambda wc: "UNDEFINED" if get_control_tag else [],
    output:
        peaks="seqnado_output/peaks/homer/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/homer/{sample}_{treatment}.log",
    params:
        options=check_options(config["homer"]["findpeaks"]),
    threads: 1
    resources:
        mem="4GB",
        runtime="2h",
    shell:
        """
        findPeaks {input.treatment} {params.options} -o {output.peaks}.tmp > {log} 2>&1 &&
        pos2bed.pl {output.peaks}.tmp -o {output.peaks} >> {log} 2>&1 &&
        rm {output.peaks}.tmp
        """


rule lanceotron_with_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=get_control_bigwig,
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.log",
    params:
        threshold=get_lanceotron_threshold,
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    threads: 1
    resources:
        mem="10GB",
        runtime="6h",
    shell:
        """
        lanceotron callPeaksInput {input.treatment} -i {input.control} -f {params.outdir} --skipheader > {log} 2>&1 &&
        cat {params.basename}_L-tron.bed | awk 'BEGIN{{OFS="\\t"}} $4 >= {params.threshold} {{print $1, $2, $3}}' > {output.peaks} 
        """


rule lanceotron_no_input:
    input:
        treatment="seqnado_output/bigwigs/deeptools/unscaled/{sample}_{treatment}.bigWig",
        control=lambda wc: "UNDEFINED" if get_control_bigwig else [],
    output:
        peaks="seqnado_output/peaks/lanceotron/{sample}_{treatment}.bed",
    log:
        "seqnado_output/logs/lanceotron/{sample}_{treatment}.bed",
    params:
        options=check_options(config["lanceotron"]["callpeak"]),
        outdir=lambda wc, output: os.path.dirname(output.peaks),
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
    container:
        "library://asmith151/seqnado/seqnado_extra:latest"
    resources:
        mem=10_1000,
        runtime="6h",
    shell:
        """
        lanceotron callPeaks {input.treatment} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
        cat {params.basename}_L-tron.bed | cut -f 1-3 > {output.peaks}
        """


ruleorder: lanceotron_with_input > lanceotron_no_input > homer_with_input > homer_no_input > macs2_with_input > macs2_no_input
