#!/usr/bin/env python3

import sys
from SigProfilerSimulator.controllers import cli_controller


def main_function():
    commands = {
        "simulate_mutations": "Simulate mutations from a mutational catalogue using genome reference files.",
    }

    if len(sys.argv) < 2 or sys.argv[1].lower() not in commands:
        _print_usage(commands)
        return

    command = sys.argv[1].lower()
    args = sys.argv[2:]

    controller = cli_controller.CliController()

    if command == "simulate_mutations":
        controller.dispatch_mutational_simulation(args)


def _print_usage(commands):
    print("Usage: SigProfilerSimulator <command> [<args>]\n")
    print("Commands:")
    for cmd, desc in commands.items():
        print(f"  {cmd}: {desc}")


if __name__ == "__main__":
    main_function()
