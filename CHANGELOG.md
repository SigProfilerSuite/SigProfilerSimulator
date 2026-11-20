# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.0] - 2026-11-20

### Added

- feat: add CLI manager with simulate_mutations command


### Changed

- docs: add MkDocs documentation site

- Documentation updated

- Cache files have been removed

- Documentation link added and dcumentacion badge corrected

- Version has been corrected

- compilation test has been activated

- tests has been disabled before create them

- bug: DBS-78 context check has been corrected

- test: First test suite has been included

- test: add unit and file-based tests for mutational_simulator pure functions

- test: integration tests has been implemented

- test: test coverage improved

- refactor: replace circular self-import with direct version import

- refactor: remove unused imports argparse, subprocess, datetime

- test: add reproducibility tests for stochastic functions

- test: add reproducibility tests for full simulation with fixed seed

- refactor: remove duplicate reference_paths() calls and dead os.system comments

- test: expand integration tests for noise, multi-sim, gender, seqInfo, DBS186, ID415

- style: convert tabs to 4-space indentation in SigProfilerSimulator.py

- style: remove dead profiling import comments from mutational_simulator.py


### Fixed

- fix: correct exome bed_file path and ref_dir shadowing bug in SigProfilerSimulator

- fix: use fork context instead of spawn to avoid forkserver crash and slowdown

- fix: probability_mask fails for chromosomes with a single mask region

- fix(simulator): fix sys.stderr.write() TypeError in INDEL exception handlers

- fix: remove module-level timer, fix NameError on stalled runs, narrow bare except

- fix: replace sys.exit() with ValueError in library functions

- fix: remove phantom pool and track combine_simulation_files results

- fix: remove dead commented block and bare pass in probability()

- fix: make seed_file parameter deterministic; default behavior unchanged

- fix: move if updating: blocks before break in INDEL simulator loops


### Fixed

- fix:multiprocessing mode changed to spawn

# Changelog
## [1.2.3] - 2026-05-06

### 💼 Other

- README has been reduced
- Automatic release has been added
- Changelogs are now aoutomated
- Merge pull request #16 from SigProfilerSuite/dev

Readme has ben reduced and release has been automated
- Revert "Merge pull request #13 from AlexandrovLab/seeds_support"

This reverts commit 596edd02d50fc37dd13505da99486a1abfd14a46, reversing
changes made to 1e8bb243d24abf20817ca8e87e64d2577a71b019.
- Revert "Merge pull request #13 from AlexandrovLab/seeds_support"

This reverts commit 596edd02d50fc37dd13505da99486a1abfd14a46, reversing
changes made to 1e8bb243d24abf20817ca8e87e64d2577a71b019.
- Merge branch 'master' into dev
- GH workflow updated
- Merge pull request #17 from SigProfilerSuite/dev

Readme has ben reduced and release has been automated
- Version updated
- Merge pull request #18 from SigProfilerSuite/dev

version updated

### 🧪 Testing

- Tests changed to GHA


All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.0] - 2025-02-24

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**