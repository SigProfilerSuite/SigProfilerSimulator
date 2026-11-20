import argparse
from typing import List

from SigProfilerSimulator.SigProfilerSimulator import SigProfilerSimulator


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def str2list(arg):
    return arg.split(",")


def parse_arguments_simulate_mutations(args: List[str], description: str) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=description)

    # Positional (required)
    parser.add_argument("project", help="Unique name for the project.")
    parser.add_argument("project_path", help="Path to the input catalogue files.")
    parser.add_argument("genome", help="Reference genome (e.g. GRCh37, GRCh38, mm10).")
    parser.add_argument("contexts", type=str2list, help="Mutation contexts, comma-separated (e.g. SBS96,DBS78,ID83).")

    # Optional
    parser.add_argument("--exome", type=str2bool, nargs="?", const=True, default=False, help="Restrict simulation to exonic regions (default: False).")
    parser.add_argument("--simulations", type=int, default=1, help="Number of simulations per sample (default: 1).")
    parser.add_argument("--updating", type=str2bool, nargs="?", const=True, default=False, help="Update chromosome sequence after each mutation (default: False).")
    parser.add_argument("--bed_file", default=None, help="BED file with genomic regions to restrict simulation.")
    parser.add_argument("--overlap", type=str2bool, nargs="?", const=True, default=False, help="Allow overlapping mutations (default: False).")
    parser.add_argument("--gender", default="female", choices=["female", "male"], help="Sample gender (default: female).")
    parser.add_argument("--seqInfo", type=str2bool, nargs="?", const=True, default=False, help="Output sequence context information (default: False).")
    parser.add_argument("--chrom_based", type=str2bool, nargs="?", const=True, default=False, help="Run simulation per chromosome (default: False).")
    parser.add_argument("--seed_file", default=None, help="Path to a seeds file for reproducible runs.")
    parser.add_argument("--spacing", type=int, default=1, help="Minimum spacing between mutations (default: 1).")
    parser.add_argument("--noisePoisson", type=str2bool, nargs="?", const=True, default=False, help="Add Poisson noise to mutation counts (default: False).")
    parser.add_argument("--noiseUniform", type=int, default=0, help="Add uniform noise (percentage) to mutation counts (default: 0).")
    parser.add_argument("--cushion", type=int, default=100, help="Cushion around chromosome ends (default: 100).")
    parser.add_argument("--region", default=None, help="Restrict simulation to a single chromosome (e.g. '1').")
    parser.add_argument("--vcf", type=str2bool, nargs="?", const=True, default=False, help="Output VCF files per sample (default: False).")
    parser.add_argument("--mask", default=None, help="Path to a mutation rate mask file.")

    return parser.parse_args(args)


class CliController:
    def dispatch_mutational_simulation(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_simulate_mutations(
            user_args, "Simulate mutations from a mutational catalogue using genome reference files."
        )
        SigProfilerSimulator(
            project=parsed_args.project,
            project_path=parsed_args.project_path,
            genome=parsed_args.genome,
            contexts=parsed_args.contexts,
            exome=parsed_args.exome,
            simulations=parsed_args.simulations,
            updating=parsed_args.updating,
            bed_file=parsed_args.bed_file,
            overlap=parsed_args.overlap,
            gender=parsed_args.gender,
            seqInfo=parsed_args.seqInfo,
            chrom_based=parsed_args.chrom_based,
            seed_file=parsed_args.seed_file,
            spacing=parsed_args.spacing,
            noisePoisson=parsed_args.noisePoisson,
            noiseUniform=parsed_args.noiseUniform,
            cushion=parsed_args.cushion,
            region=parsed_args.region,
            vcf=parsed_args.vcf,
            mask=parsed_args.mask,
        )
