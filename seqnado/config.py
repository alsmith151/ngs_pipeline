import datetime
import json
import os
import pathlib
import sys

from jinja2 import Environment, FileSystemLoader
from loguru import logger

logger.add(sys.stderr, level="INFO")

package_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.join(package_dir, "workflow/config")


def get_user_input(prompt, default=None, is_boolean=False, choices=None):
    while True:
        user_input = (
            input(f"{prompt} [{'/'.join(choices) if choices else default}]: ")
            or default
        )
        if is_boolean:
            return user_input.lower() == "yes"
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from {', '.join(choices)}.")
            continue
        return user_input


def setup_configuration(assay, template_data, seqnado_version):
    genome_config_path = pathlib.Path(
        os.getenv(
            "SEQNADO_CONFIG",
            pathlib.Path.home(),
        )
    )
    genome_config_file = (
        genome_config_path / ".config" / "seqnado" / "genome_config.json"
    )
    if not genome_config_file.exists():
        logger.info(
            "Genome config file not found. Please run 'seqnado-init' to create the genome config file."
        )
        sys.exit(1)
    else:
        with open(genome_config_file) as f:
            genome_values = json.load(f)
    username = os.getenv("USER", "unknown_user")
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    project_name = get_user_input(
        "What is your project name?", default=f"{username}_project"
    ).replace(" ", "_")
    genome = get_user_input("What is the genome?", default="hg38")
    if genome in genome_values:
        genome_config = {
            "index": genome_values[genome][
                "star_index" if assay == "rna" else "bt2_index"
            ],
            "chromosome_sizes": genome_values[genome]["chromosome_sizes"],
            "gtf": genome_values[genome]["gtf"],
            "blacklist": genome_values[genome].get("blacklist", ""),
        }
        template_data.update(genome_config)
    else:
        logger.error(
            f"Genome '{genome}' not found in genome config file. Please update the genome config file: {genome_config_file}"
        )
        sys.exit(1)
    common_config = {
        "username": username,
        "project_date": today,
        "project_name": project_name,
        "seqnado_version": seqnado_version,
        "genome": genome,
    }
    template_data.update(common_config)
    # Fastqscreen
    template_data["fastq_screen"] = get_user_input(
        "Perform fastqscreen? (yes/no)", default="no", is_boolean=True
    )
    if template_data["fastq_screen"]:
        template_data["fastq_screen_config"] = get_user_input(
            "Path to fastqscreen config:",
            default="/ceph/project/milne_group/shared/seqnado_reference/fastqscreen_reference/fastq_screen.conf",
        )
    # Blacklist
    template_data["remove_blacklist"] = get_user_input(
        "Do you want to remove blacklist regions? (yes/no)",
        default="yes",
        is_boolean=True,
    )
    if template_data["remove_blacklist"]:
        template_data["blacklist"] = genome_values[genome]["blacklist"]
    # Handle duplicates
    template_data["remove_pcr_duplicates"] = get_user_input(
        "Remove PCR duplicates? (yes/no)",
        default="yes" if assay in ["chip", "atac"] else "no",
        is_boolean=True,
    )
    if template_data["remove_pcr_duplicates"]:
        template_data["remove_pcr_duplicates_method"] = get_user_input(
            "Remove PCR duplicates method:",
            default="picard",
            choices=["picard", "samtools"],
        )
        # Library Complexity
        template_data["library_complexity"] = get_user_input(
            "Calculate library complexity? (yes/no)", default="no", is_boolean=True
        )
    else:
        template_data["remove_pcr_duplicates_method"] = "False"
        template_data["library_complexity"] = "False"
    # Shift reads
    if assay == "atac":
        template_data["shift_atac_reads"] = (
            get_user_input(
                "Shift ATAC-seq reads? (yes/no)", default="yes", is_boolean=True
            )
            if assay == "atac"
            else "False"
        )
    # Spike in
    if assay in ["chip", "rna"]:
        template_data["spikein"] = get_user_input(
            "Do you have spikein? (yes/no)", default="no", is_boolean=True
        )
        if template_data["spikein"] and not assay == "rna":
            template_data["normalisation_method"] = get_user_input(
                "Normalisation method:",
                default="orlando",
                choices=["orlando", "with_input"],
            )
            template_data["reference_genome"] = get_user_input(
                "Reference genome:", default="hg38"
            )
            template_data["spikein_genome"] = get_user_input(
                "Spikein genome:", default="dm6"
            )
    # Make bigwigs
    if assay not in ["snp"]:
        template_data["make_bigwigs"] = get_user_input(
            "Do you want to make bigwigs? (yes/no)", default="no", is_boolean=True
        )
        if template_data["make_bigwigs"]:
            template_data["pileup_method"] = get_user_input(
                "Pileup method:",
                default="deeptools",
                choices=["deeptools", "homer"],
            )
            # template_data["scale"] = get_user_input(
            #     "Scale bigwigs? (yes/no)", default="no", is_boolean=True
            # )
            template_data["make_heatmaps"] = get_user_input(
                "Do you want to make heatmaps? (yes/no)", default="no", is_boolean=True
            )
        else:
            template_data["pileup_method"] = "False"
            template_data["scale"] = "False"
            template_data["make_heatmaps"] = "False"
    # Call peaks
    if assay in ["chip", "atac"]:
        template_data["call_peaks"] = get_user_input(
            "Do you want to call peaks? (yes/no)", default="no", is_boolean=True
        )
        if template_data["call_peaks"]:
            template_data["peak_calling_method"] = get_user_input(
                "Peak caller:",
                default="lanceotron",
                choices=["lanceotron", "macs", "homer", "seacr"],
            )
    if assay in ["chip", "atac"]:
        template_data["consenus_counts"] = get_user_input(
            "Generate consensus counts from Design merge column? (yes/no)",
            default="no",
            is_boolean=True,
        )
    # RNA options
    template_data["rna_quantification"] = (
        get_user_input(
            "RNA quantification method:",
            default="feature_counts",
            choices=["feature_counts", "salmon"],
        )
        if assay == "rna"
        else "False"
    )
    template_data["salmon_index"] = (
        get_user_input(
            "Path to salmon index:",
            default="path/to/salmon_index",
        )
        if template_data["rna_quantification"] == "salmon"
        else "False"
    )
    # Run DESeq2
    template_data["run_deseq2"] = (
        get_user_input("Run DESeq2? (yes/no)", default="no", is_boolean=True)
        if assay == "rna"
        else "False"
    )
    # SNP options
    template_data["call_snps"] = (
        get_user_input("Call SNPs? (yes/no)", default="no", is_boolean=True)
        if assay == "snp"
        else "False"
    )
    if assay == "snp" and template_data["call_snps"]:
        template_data["snp_calling_method"] = get_user_input(
            "SNP caller:",
            default="bcftools",
            choices=["bcftools", "deepvariant"],
        )
        template_data["fasta"] = get_user_input(
            "Path to reference fasta:", default="path/to/reference.fasta"
        )
        template_data["fasta_index"] = get_user_input(
            "Path to reference fasta index:", default="path/to/reference.fasta.fai"
        )
        template_data["snp_database"] = get_user_input(
            "Path to SNP database:",
            default="path/to/snp_database",
        )
    else:
        template_data["snp_calling_method"] = "False"
        template_data["fasta"] = "False"
        template_data["fasta_index"] = "False"
        template_data["snp_database"] = "False"
    # Make UCSC hub
    template_data["make_ucsc_hub"] = get_user_input(
        "Do you want to make a UCSC hub? (yes/no)", default="no", is_boolean=True
    )
    template_data["UCSC_hub_directory"] = (
        get_user_input("UCSC hub directory:", default="seqnado_output/hub/")
        if template_data["make_ucsc_hub"]
        else "seqnado_output/hub/"
    )
    template_data["email"] = (
        get_user_input("What is your email address?", default=f"{username}@example.com")
        if template_data["make_ucsc_hub"]
        else f"{username}@example.com"
    )
    template_data["color_by"] = (
        get_user_input("Color by (for UCSC hub):", default="samplename")
        if template_data["make_ucsc_hub"]
        else "samplename"
    )
    template_data["options"] = (
        TOOL_OPTIONS
        if assay in ["chip", "atac"]
        else (
            TOOL_OPTIONS_RNA
            if assay == "rna"
            else TOOL_OPTIONS_SNP
            if assay == "snp"
            else ""
        )
    )
    template_data["geo_submission_files"] = get_user_input(
        "Generate GEO submission files (MD5Sums, read count summaries...)? (yes/no)",
        default="no",
        is_boolean=True,
    )
    template_data["perform_plotting"] = get_user_input(
        "Perform plotting? (yes/no)", default="no", is_boolean=True
    )
    if template_data["perform_plotting"]:
        template_data["plotting_coordinates"] = get_user_input(
            "Path to bed file with coordinates for plotting", default=None
        )
        if genome in genome_values and genome_values.get("genes"):
            template_data["plotting_genes"] = genome_values[genome].get("genes")
        else:
            template_data["plotting_genes"] = get_user_input(
                "Path to bed file with genes.", default=None
            )
    else:
        template_data["plotting_coordinates"] = None
        template_data["plotting_genes"] = None


TOOL_OPTIONS = """
trim_galore:
    threads: 4
    options: --2colour 20 

bowtie2:
    threads: 8
    options:

samtools:
    threads: 16
    filter_options: -f 2

picard:
    threads: 8
    options:

homer:
    use_input: true
    maketagdirectory:
    makebigwig:
    findpeaks:

deeptools:
    threads: 8
    alignmentsieve: --minMappingQuality 30 
    bamcoverage: --extendReads -bs 1 --normalizeUsing RPKM --minMappingQuality 10

macs:
    version: 2
    callpeak: -f BAMPE

lanceotron:
    use_input: True
    callpeak: -c 0.5

seacr:
    threshold: 0.01
    norm: non
    stringency: stringent

heatmap:
    options: -b 1000 -m 5000 -a 1000
    colormap: RdYlBu_r 

featurecounts:
    threads: 16
    options:  -p --countReadPairs
    
"""

TOOL_OPTIONS_RNA = """
trim_galore:
    threads: 4
    options: --2colour 20 

star:
    threads: 16
    options: --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard

samtools:
    threads: 16
    filter_options: -f 2

picard:
    threads: 8
    options:

featurecounts:
    threads: 16
    options: -s 0 -p --countReadPairs -t exon -g gene_id

salmon:
    threads: 16
    options: --libType A
    
homer:
    maketagdirectory:
    makebigwig:

deeptools:
    threads: 16
    alignmentsieve: --minMappingQuality 30 
    bamcoverage: -bs 1 --normalizeUsing CPM

heatmap:
    options: -b 1000 -m 5000 -a 1000
    colormap: RdYlBu_r 
"""

TOOL_OPTIONS_SNP = """
trim_galore:
    threads: 8
    options: --2colour 20 

bowtie2:
    threads: 8
    options:

samtools:
    threads: 16
    filter_options: -f 2
    
picard:
    threads: 8
    options:

bcftools:
    threads: 16
    options:
    
"""


def create_config(assay, rerun, seqnado_version, debug=False):
    env = Environment(loader=FileSystemLoader(template_dir), auto_reload=False)
    template = env.get_template("config.yaml.jinja")
    template_deseq2 = env.get_template("deseq2.qmd.jinja")
    # Initialize template data
    template_data = {"assay": assay, "seqnado_version": seqnado_version}
    # Setup configuration
    setup_configuration(assay, template_data, seqnado_version)
    # Create directory and render template
    if rerun:
        dir_name = os.getcwd()
        with open(os.path.join(dir_name, f"config_{assay}.yml"), "w") as file:
            file.write(template.render(template_data))
    else:
        dir_name = f"{template_data['project_date']}_{template_data['assay']}_{template_data['project_name']}"
        os.makedirs(dir_name, exist_ok=True)
        fastq_dir = os.path.join(dir_name, "fastq")
        os.makedirs(fastq_dir, exist_ok=True)
        with open(os.path.join(dir_name, f"config_{assay}.yml"), "w") as file:
            file.write(template.render(template_data))
    # add deseq2 qmd file if rna
    if assay == "rna":
        with open(
            os.path.join(dir_name, f"deseq2_{template_data['project_name']}.qmd"), "w"
        ) as file:
            file.write(template_deseq2.render(template_data))
    print(
        f"Directory '{dir_name}' has been created with the 'config_{assay}.yml' file."
    )
    if debug:
        with open(os.path.join(dir_name, "data.json"), "w") as file:
            json.dump(template_data, file)
