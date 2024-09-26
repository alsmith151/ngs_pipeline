import click
import os
import subprocess
import re
from loguru import logger
import sys
import pathlib
import shlex


FILE = os.path.abspath(__file__)
PACKAGE_DIR = os.path.dirname(FILE)


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument("method", type=click.Choice(["atac", "chip",'cat',"rna", "snp"]))
@click.option("-r", "--rerun", is_flag=True, help="Re-run the config")
@click.option(
    "-g",
    "--genome",
    default="other",
    help="Genome to use",
    type=click.Choice(
        choices=[
            "dm6",
            "hg19",
            "hg38",
            "hg38_dm6",
            "hg38_mm39",
            "hg38_spikein",
            "mm10",
            "mm39",
            "other",
        ]
    ),
)
def cli_config(method, help=False, genome="other", rerun=False):
    """
    Runs the config for the data processing pipeline.
    """
    import seqnado.config as config

    config.create_config(method, genome, rerun)


@click.command()
@click.argument("method", type=click.Choice(["atac", "chip", "rna",'cat',"snp"]))
@click.argument("files", nargs=-1)
@click.option("-o", "--output", default="design.csv", help="Output file name")
def cli_design(method, files, output="design.csv"):
    """
    Generates a SeqNado design file from a list of files.
    """
    import pathlib
    from seqnado.design import Design, DesignIP, FastqFile, FastqFileIP
    
    if not files:
        potential_file_locations = [
            ".",
            "fastqs",
            "fastq",
            "data",
            "data/fastqs",
        ]

        for location in potential_file_locations:
            files = list(pathlib.Path(location).glob("*.fastq.gz"))
            if files:
                break

        if not files:
            logger.error("No fastq files provided or found in current directory")
            logger.error(f"""
                         Fastq files can be provided as arguments or found in the following directories:
                         {potential_file_locations}
                         """)
            raise ValueError("No fastq files provided or found in current directory" )



    if not method == "chip":
        design = Design.from_fastq_files(files)
    else:

        design = DesignIP.from_fastq_files(files)

    (
        design.to_dataframe()
        .assign(scale_group="all")
        .sort_values("sample_name")
        .to_csv(output, index=False)
    )


@click.command(context_settings=dict(ignore_unknown_options=True))
@click.argument(
    "method",
    type=click.Choice(["atac", "chip", 'cat', "rna", "snp"]),
)
@click.option("--version", help="Print version and exit", is_flag=True)
@click.option(
    "--preset",
    default="lc",
    help="""Pre-set snakemake job profile to use for pipeline run:
            lc: local conda environment
            ls: local singularity environment
            ss: slurm singularity environment (runs jobs on cluster)
            """,
    type=click.Choice(choices=["lc", "ls", "ss"]),
)
@click.option(
    "--clean-symlinks",
    is_flag=True,
    help="Remove symlinks created by previous runs. Useful for re-running pipeline after misconfiguration.",
)
@click.option(
    '-s',
    '--scale-resources',
    help="Scale factor the memory and time resources for the pipeline",
    default=1.0,
    type=float
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Increase logging verbosity",
)
@click.argument("pipeline_options", nargs=-1, type=click.UNPROCESSED)
def cli_pipeline(
    method,
    pipeline_options,
    help=False,
    preset="local",
    version=False,
    verbose=False,
    clean_symlinks=False,
    scale_resources=1.0,
):
    """Runs the data processing pipeline"""

    from seqnado.helpers import extract_cores_from_options

    if version:
        from importlib.metadata import version

        _version = version("seqnado")

        _version = version("seqnado")
        print(f"SeqNado version {_version}")
        sys.exit(0)
    
    if verbose:
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
    else:
        logger.remove()
        logger.add(sys.stderr, level="INFO")

    pipeline_options, cores = extract_cores_from_options(pipeline_options)

    # Scale the memory and time resources
    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    # Removes old symlinks if requested
    if clean_symlinks:
        logger.info("Cleaning symlinks")
        links = pathlib.Path("seqnado_output/fastqs").glob("*")
        for link in links:
            if link.is_symlink():
                link.unlink()
    
    cmd = [
        "snakemake",
        "-c",
        str(cores),
        "--snakefile",
        os.path.join(PACKAGE_DIR, "workflow", f"snakefile_{method.replace('-', '_')}"),
    ]

    if pipeline_options:
        cmd.extend(pipeline_options)

    if preset == "ss":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(
                    os.path.join(
                        PACKAGE_DIR, "workflow/envs/profiles/profile_slurm_singularity"
                    )
                ),
            ]
        )

    elif preset == "ls":
        cmd.extend(
            [
                "--profile",
                os.path.abspath(
                    os.path.join(
                        PACKAGE_DIR, "workflow/envs/profiles/profile_singularity"
                    )
                ),
            ]
        )

    cmd.extend(["--show-failed-logs"])

    # Print the logo
    with open(f"{PACKAGE_DIR}/data/logo.txt", "r") as f:
        logo = f.read()

    print(logo)



    # Home directory symlinks cause issues with singularity bind mounts
    # to avoid this will change directory to the full resolved path of the current directory
    cwd = str(pathlib.Path(".").resolve())
    os.chdir(cwd)
    os.environ["PWD"] = cwd
    completed = subprocess.run(cmd, cwd=cwd)

    if completed.returncode != 0:
        sys.exit(completed.returncode)
    else:
        sys.exit(0)
