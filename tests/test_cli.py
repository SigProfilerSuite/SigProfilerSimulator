"""
Unit tests for the CLI controller.
Verifies that CLI arguments are parsed and forwarded to SigProfilerSimulator correctly.
No genome reference files required.
Run with: pytest tests/test_cli.py -v
"""

import sys
import os
from unittest.mock import patch

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from SigProfilerSimulator.controllers.cli_controller import CliController


_MOCK_TARGET = "SigProfilerSimulator.controllers.cli_controller.SigProfilerSimulator"

# Base positional args shared across tests
_BASE = ["myproj", "/data/input", "GRCh37", "SBS96"]


class TestDispatchMutationalSimulation:
    def test_required_args_forwarded(self):
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(_BASE)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["project"] == "myproj"
        assert kwargs["project_path"] == "/data/input"
        assert kwargs["genome"] == "GRCh37"
        assert kwargs["contexts"] == ["SBS96"]

    def test_multiple_contexts_comma_separated(self):
        args = ["p", "/x", "GRCh38", "SBS96,DBS78,ID83"]
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(args)
        assert mock_sim.call_args.kwargs["contexts"] == ["SBS96", "DBS78", "ID83"]

    def test_defaults(self):
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(_BASE)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["exome"] is False
        assert kwargs["simulations"] == 1
        assert kwargs["updating"] is False
        assert kwargs["bed_file"] is None
        assert kwargs["overlap"] is False
        assert kwargs["gender"] == "female"
        assert kwargs["seqInfo"] is False
        assert kwargs["chrom_based"] is False
        assert kwargs["seed_file"] is None
        assert kwargs["spacing"] == 1
        assert kwargs["noisePoisson"] is False
        assert kwargs["noiseUniform"] == 0
        assert kwargs["cushion"] == 100
        assert kwargs["region"] is None
        assert kwargs["vcf"] is False
        assert kwargs["mask"] is None

    def test_bool_flags_with_true_value(self):
        args = _BASE + [
            "--exome", "true", "--updating", "true", "--overlap", "true",
            "--seqInfo", "true", "--chrom_based", "true",
            "--noisePoisson", "true", "--vcf", "true",
        ]
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(args)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["exome"] is True
        assert kwargs["updating"] is True
        assert kwargs["overlap"] is True
        assert kwargs["seqInfo"] is True
        assert kwargs["chrom_based"] is True
        assert kwargs["noisePoisson"] is True
        assert kwargs["vcf"] is True

    def test_bool_flags_bare_without_value(self):
        args = _BASE + ["--exome", "--noisePoisson", "--vcf"]
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(args)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["exome"] is True
        assert kwargs["noisePoisson"] is True
        assert kwargs["vcf"] is True

    def test_bool_flags_with_false_value(self):
        args = _BASE + ["--exome", "false", "--noisePoisson", "0"]
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(args)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["exome"] is False
        assert kwargs["noisePoisson"] is False

    def test_optional_values(self):
        args = _BASE + [
            "--simulations", "5",
            "--gender", "male",
            "--spacing", "10",
            "--noiseUniform", "20",
            "--cushion", "50",
            "--region", "1",
            "--bed_file", "/data/regions.bed",
            "--seed_file", "/data/seeds.txt",
            "--mask", "/data/mask.csv",
        ]
        with patch(_MOCK_TARGET) as mock_sim:
            CliController().dispatch_mutational_simulation(args)
        kwargs = mock_sim.call_args.kwargs
        assert kwargs["simulations"] == 5
        assert kwargs["gender"] == "male"
        assert kwargs["spacing"] == 10
        assert kwargs["noiseUniform"] == 20
        assert kwargs["cushion"] == 50
        assert kwargs["region"] == "1"
        assert kwargs["bed_file"] == "/data/regions.bed"
        assert kwargs["seed_file"] == "/data/seeds.txt"
        assert kwargs["mask"] == "/data/mask.csv"

    def test_missing_required_positional_exits(self):
        with patch(_MOCK_TARGET):
            with pytest.raises(SystemExit):
                CliController().dispatch_mutational_simulation(["myproj", "/x", "GRCh37"])

    def test_invalid_gender_exits(self):
        args = _BASE + ["--gender", "other"]
        with patch(_MOCK_TARGET):
            with pytest.raises(SystemExit):
                CliController().dispatch_mutational_simulation(args)

    def test_invalid_bool_value_exits(self):
        args = _BASE + ["--exome", "maybe"]
        with patch(_MOCK_TARGET):
            with pytest.raises(SystemExit):
                CliController().dispatch_mutational_simulation(args)


class TestMainFunction:
    def test_no_args_prints_usage(self, capsys):
        with patch("sys.argv", ["SigProfilerSimulator"]):
            from SigProfilerSimulator.sigprofilesimulator_cli import main_function
            main_function()
        assert "simulate_mutations" in capsys.readouterr().out

    def test_unknown_command_prints_usage(self, capsys):
        with patch("sys.argv", ["SigProfilerSimulator", "unknown_cmd"]):
            from SigProfilerSimulator.sigprofilesimulator_cli import main_function
            main_function()
        assert "simulate_mutations" in capsys.readouterr().out

    def test_simulate_mutations_dispatches(self):
        argv = ["SigProfilerSimulator", "simulate_mutations", "p", "/x", "GRCh37", "SBS96"]
        with patch("sys.argv", argv), patch(_MOCK_TARGET) as mock_sim:
            from SigProfilerSimulator.sigprofilesimulator_cli import main_function
            main_function()
        mock_sim.assert_called_once()
