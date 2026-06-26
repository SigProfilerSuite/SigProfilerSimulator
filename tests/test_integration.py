"""
Integration tests for SigProfilerSimulator.
Require GRCh37 genome installed via SigProfilerMatrixGenerator.
Skipped automatically when the genome is not available.
"""

import os
import shutil
import tempfile

import pytest

from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as matRef
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator.SigProfilerSimulator import SigProfilerSimulator

PROJECT = "test_sps"
GENOME = "GRCh37"
CONTEXTS = ["96"]
SIMULATIONS = 1
N_INPUT_MUTATIONS = 12

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
VCF_PATH = os.path.join(TESTS_DIR, "data", "GRCh37_test.vcf")

MAF_HEADER = "\t".join([
    "Hugo_symbol", "Entrez_gene_ID", "Center", "Genome", "Chrom",
    "Start_position", "End_position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "matGenClass",
])


def _grch37_available():
    try:
        path, _ = matRef.reference_paths("GRCh37")
        return os.path.exists(path) and len(os.listdir(path)) > 2
    except Exception:
        return False


grch37 = pytest.mark.skipif(
    not _grch37_available(),
    reason="GRCh37 genome not installed",
)


@pytest.fixture(scope="module")
def simulation_output():
    """Runs a full simulation once and returns the output MAF path."""
    project_path = tempfile.mkdtemp() + "/"
    try:
        input_dir = project_path + "input/"
        os.makedirs(input_dir)
        shutil.copy(VCF_PATH, input_dir)

        matGen.SigProfilerMatrixGeneratorFunc(PROJECT, GENOME, project_path, plot=False)
        SigProfilerSimulator(PROJECT, project_path, GENOME, CONTEXTS, simulations=SIMULATIONS)

        context_str = "_".join(CONTEXTS)
        maf = (
            f"{project_path}output/simulations/"
            f"{PROJECT}_simulations_{GENOME}_{context_str}/{SIMULATIONS}.maf"
        )
        yield maf
    finally:
        shutil.rmtree(project_path, ignore_errors=True)


@grch37
class TestEndToEndSimulation:
    def test_output_maf_created(self, simulation_output):
        assert os.path.exists(simulation_output)

    def test_output_maf_header(self, simulation_output):
        with open(simulation_output) as f:
            header = f.readline().strip()
        assert header == MAF_HEADER

    def test_output_maf_mutation_count(self, simulation_output):
        with open(simulation_output) as f:
            n = sum(1 for _ in f) - 1  # subtract header
        assert n == N_INPUT_MUTATIONS

    def test_output_maf_genome_field(self, simulation_output):
        with open(simulation_output) as f:
            f.readline()
            first_data = f.readline().strip().split("\t")
        assert first_data[3] == GENOME

    def test_output_maf_variant_type_snp(self, simulation_output):
        with open(simulation_output) as f:
            lines = f.readlines()[1:]
        variant_types = {line.strip().split("\t")[9] for line in lines if line.strip()}
        assert variant_types == {"SNP"}
