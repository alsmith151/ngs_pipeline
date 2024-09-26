import glob
import os
import pathlib
import re
import shutil
import subprocess
import tarfile
from datetime import datetime
import shlex

import pytest
import requests


@pytest.fixture(
    scope="function",
    params=["atac", "chip", "chip-rx", "rna", "rna-rx", "snp", "cat"],
    autouse=True,
)
def assay(request):
    return request.param


@pytest.fixture(scope="function")
def repo_path() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="function")
def package_path(repo_path):
    return repo_path / "seqnado"


@pytest.fixture(scope="function")
def test_dir_path(repo_path):
    return repo_path / "tests"


@pytest.fixture(scope="function")
def test_data_path(test_dir_path):
    return test_dir_path / "data"


@pytest.fixture(scope="function")
def workflow_path(package_path):
    return package_path / "workflow"


@pytest.fixture(scope="function")
def config_path(workflow_path):
    return workflow_path / "config"


@pytest.fixture(scope="function")
def genome_path(test_data_path):
    p = test_data_path / "genome"
    p.mkdir(parents=True, exist_ok=True)
    return p


@pytest.fixture(scope="function")
def genome_indices_path(genome_path, assay) -> pathlib.Path:

    if "rna" not in assay:
        return genome_path / "bt2_chr21_dm6_chr2L"
    else:
        return genome_path / "STAR_chr21_rna_spikein"


@pytest.fixture(scope="function")
def indicies(genome_indices_path, genome_path) -> pathlib.Path:

    download_indices = True if not genome_indices_path.exists() else False
    suffix = genome_indices_path.with_suffix(".tar.gz").name
    url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/{suffix}"

    if "bt2" in str(genome_indices_path):
        indicies_path = genome_indices_path / "bt2_chr21_dm6_chr2L"
    else:
        indicies_path = genome_indices_path

    if download_indices:

        r = requests.get(url, stream=True)

        tar_index = genome_indices_path.with_suffix(".tar.gz")

        with open(tar_index, "wb") as f:
            f.write(r.content)

        tar = tarfile.open(tar_index)

        if "bt2" in str(
            genome_indices_path
        ):  # These are individual files so need to extract to the indicies folder
            genome_indices_path.mkdir(parents=True, exist_ok=True)
            tar.extractall(path=genome_indices_path, filter="data")
        else:
            tar.extractall(path=genome_path, filter="data")
        tar.close()

        os.remove(tar_index)

    return indicies_path


@pytest.fixture(scope="function")
def chromsizes(genome_path):

    suffix = "chr21.fa.fai"

    if not (genome_path / "chr21_rename.fa.fai").exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{suffix}"
        r = requests.get(url, stream=True)
        with open(genome_path / suffix, "wb") as f:
            f.write(r.content)

    return genome_path / suffix


@pytest.fixture(scope="function")
def gtf(genome_path, assay, indicies):

    if "rna" in assay:
        gtf_path = genome_path / "chr21_rna_spikein.gtf"
    else:
        gtf_path = genome_path / "chr21.gtf"

    if not gtf_path.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/asmith/ngs_pipeline/{gtf_path.name}"
        r = requests.get(url, stream=True)
        with open(gtf_path, "wb") as f:
            f.write(r.content)

    return gtf_path


@pytest.fixture(scope="function")
def blacklist(genome_path):

    blacklist_path = genome_path / "hg38-blacklist.v2.bed.gz"
    if not blacklist_path.exists():
        url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
        r = requests.get(url, stream=True)
        with open(blacklist_path, "wb") as f:
            f.write(r.content)

    return blacklist_path


@pytest.fixture(scope="function")
def assay_type(assay):
    return re.sub(r"(.*)\-.*", r"\1", assay)


@pytest.fixture(scope="function")
def snakefile_path(package_path, assay_type):
    return package_path / "workflow" / f"snakefile_{assay_type}"


@pytest.fixture(scope="function")
def fastqs(test_data_path, assay) -> list[pathlib.Path]:
    path = test_data_path / "fastq"

    if not path.exists():
        url = f"https://userweb.molbiol.ox.ac.uk/public/project/milne_group/cchahrou/seqnado_reference/fastq.tar.gz"
        r = requests.get(url, stream=True)

        tar_path = path.with_suffix(".tar.gz")

        with open(tar_path, "wb") as f:
            f.write(r.content)

        with tarfile.open(tar_path) as tar:
            tar.extractall(path=path.parent, filter="data")

    match assay:
        case "atac":
            files = list(path.glob("atac*.fastq.gz"))
        case "chip":
            files = list(path.glob("chip-rx*.fastq.gz"))
            files.append(path / "chip-rx-single_MLL.fastq.gz")
        case "chip-rx":
            files = list(path.glob("chip-rx*.fastq.gz"))
        case "rna":
            files = list(path.glob("rna_*.fastq.gz"))
        case "rna-rx":
            files = list(path.glob("rna-spikein*.fastq.gz"))
        case "snp":
            files = list(path.glob("snp*.fastq.gz"))

    return files


@pytest.fixture(scope="function")
def run_directory(tmpdir_factory, assay):
    fn = tmpdir_factory.mktemp(assay)
    return fn


@pytest.fixture(scope="function")
def user_inputs(
    test_data_path, indicies, chromsizes, assay, assay_type, gtf, blacklist
):

    defaults = {
        "project_name": "test",
        "genome_name": "hg38",
        "indices": str(indicies),
        "chromsizes": str(chromsizes),
        "gtf": str(gtf),
        "blacklist": str(blacklist),
        "fastq_screen": "no",
        "remove_blacklist": "yes",
    }

    defaults_atac = {
        "remove_pcr_duplicates": "yes",
        "remove_pcr_duplicates_method": "picard",
        "library_complexity": "yes",
        "shift_atac_reads": "yes",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "yes",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_chip = {
        "remove_pcr_duplicates": "no",
        "spikein": "no",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "yes",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_chip_rx = {
        "remove_pcr_duplicates": "no",
        "spikein": "yes",
        "normalisation_method": "orlando",
        "reference_genome": "hg38",
        "spikein_genome": "dm6",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "no",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
        "call_peaks": "yes",
        "peak_calling_method": "lanceotron",
    }

    defaults_cat = defaults_chip.copy()
    defaults_cat['peak_calling_method'] = 'seacr'


    defaults_rna = {
        "remove_pcr_duplicates": "no",
        "spikein": "no",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "yes",
        "rna_quantification": "feature_counts",
        "run_deseq2": "no",
    }

    defaults_rna_rx = {
        "remove_pcr_duplicates": "no",
        "spikein": "yes",
        "make_bigwigs": "yes",
        "pileup_method": "deeptools",
        "scale": "no",
        "make_heatmaps": "no",
        "rna_quantification": "feature_counts",
        "run_deseq2": "yes",
    }

    defaults_snp = {
        "remove_pcr_duplicates": "no",
        "call_snps": "no",
    }

    hub = {
        "make_ucsc_hub": "yes",
        "UCSC_hub_directory": "test_hub",
        "email": "test",
        "color_by": "samplename",
    }

    match assay:
        case "atac":
            return {**defaults, **defaults_atac, **hub}
        case "chip":
            return {**defaults, **defaults_chip, **hub}
        case "chip-rx":
            return {**defaults, **defaults_chip_rx, **hub}
        case "rna":
            return {**defaults, **defaults_rna, **hub}
        case "rna-rx":
            return {**defaults, **defaults_rna_rx, **hub}
        case "snp":
            return {**defaults, **defaults_snp, **hub}
        case "cat":
            return {**defaults, **defaults_cat, **hub}


@pytest.fixture(scope="function")
def config_yaml(run_directory, user_inputs, assay_type):

    user_inputs = "\n".join([str(v) for v in user_inputs.values()])
    cmd = ["seqnado-config", assay_type]

    # Run the script with subprocess
    process = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=run_directory,
    )

    stdout, stderr = process.communicate(input=user_inputs)

    project_name = "test"
    date = datetime.now().strftime("%Y-%m-%d")
    config_file_path = (
        run_directory / f"{date}_{assay_type}_{project_name}/config_{assay_type}.yml"
    )
    if not config_file_path.exists():
        print("debug stdout", stdout)

    assert config_file_path.exists(), f"{assay_type} config file not created."
    return config_file_path


@pytest.fixture(scope="function")
def config_yaml_for_testing(config_yaml, assay):
    import yaml

    with open(config_yaml, "r") as f:
        config = yaml.safe_load(f)

    if assay == "chip":
        config["library_complexity"] = False
    elif assay == "chip-rx":
        config["peak_calling_method"] = "seacr"
    elif assay == "cat":
        config["peak_calling_method"] = "seacr"
    elif assay == "atac":
        config["pileup_method"] = ["deeptools", "homer"]
        config["call_peaks"] = True
        config["peak_calling_method"] = ["lanceotron", "macs", "homer"]
        config["peak_calling_method"] = ["lanceotron", "macs", "homer"]

    with open(config_yaml, "w") as f:
        yaml.dump(config, f)

    return pathlib.Path(config_yaml)


@pytest.fixture(scope="function")
def seqnado_run_dir(config_yaml_for_testing):
    return pathlib.Path(config_yaml_for_testing).parent


@pytest.fixture(scope="function")
def design(seqnado_run_dir, assay_type, assay):
    cmd = ["seqnado-design", assay_type]
    completed = subprocess.run(" ".join(cmd), shell=True, cwd=seqnado_run_dir)
    assert completed.returncode == 0

    if assay == "chip":
        # Add merge column to design file
        import pandas as pd

        df = pd.read_csv(seqnado_run_dir / "design.csv", index_col=0)
        df["merge"] = "MLL-MERGED-TOGETHER"
        df.to_csv(seqnado_run_dir / "design.csv")

    elif assay == "rna-rx":
        # Add deseq2 column to design file
        import pandas as pd

        df = pd.read_csv(seqnado_run_dir / "design.csv", index_col=0)
        df["deseq2"] = df.index.str.split("-").str[-2]
        df.to_csv(seqnado_run_dir / "design.csv")

    return seqnado_run_dir / "design.csv"


@pytest.fixture(scope="function", autouse=True)
def set_up(seqnado_run_dir, fastqs):

    cwd = pathlib.Path(os.getcwd())
    os.chdir(seqnado_run_dir)

    # Move fastqs to run directory
    for fq in fastqs:
        shutil.copy(fq, seqnado_run_dir)

    yield

    os.chdir(cwd)


def test_config(config_yaml, assay_type):
    assert os.path.exists(config_yaml), f"{assay_type} config file not created."

def test_design(design, assay_type):
    assert os.path.exists(design), f"{assay_type} design file not created."

@pytest.fixture(scope="function", autouse=True)
def apptainer_args(indicies, test_data_path):
    indicies_mount = indicies.parent if not indicies.is_dir() else indicies
    tmpdir = pathlib.Path(os.environ.get("TMPDIR", "/tmp") or "/tmp")
    wd = pathlib.Path(os.getcwd()).resolve()
    apptainer_cache_dir = pathlib.Path.home() / ".apptainer"
    
    os.environ["APPTAINER_BINDPATH"] = (
        f"{wd}:{wd},"
        f"{test_data_path}:{test_data_path},"
        f"{indicies_mount}:{indicies_mount},"
        f"{tmpdir}:{tmpdir}"
    )
    os.environ["APPTAINER_CACHEDIR"] = str(apptainer_cache_dir)


@pytest.fixture(scope="function")
def test_profile_path(workflow_path):
    return workflow_path / "envs" / "profiles" / "profile_test"


def test_pipeline(
    assay_type,
    config_yaml_for_testing,
    apptainer_args,
    cores,
    design,
    indicies,
    test_data_path,
    test_profile_path,
):

    subprocess.run(
        [
            "seqnado",
            assay_type,
            "-c",
            str(cores),
            "--configfile",
            str(config_yaml_for_testing),
            "--workflow-profile",
            test_profile_path,
            "--apptainer-prefix",
            "/tmp",
        ],
    )

    assert not os.path.exists("seqnado_error.log")
    assert os.path.exists("seqnado_output/")
    assert os.path.exists("seqnado_output/qc/full_qc_report.html")