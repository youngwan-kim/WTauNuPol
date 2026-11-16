# TauPolarization

Longitudinal τ polarization in W → τν on SKNanoAnalyzer.

## Layout

| path | content |
| --- | --- |
| `include/TauPolarization/TauPolarimeter.h`, `src/TauPolarimeter.cc` | per-τ polarimeter ω, template weight, ratio helpers, shape Fisher information; pure functions |
| `include/TauPolarization/GenTauBuilder.h`, `src/GenTauBuilder.cc` | τ from W (and t → W) out of flat GenPart columns; framework-independent |
| `include/TauPolarization/TauPolarization.h`, `src/TauPolarization.cc` | the analyzer |
| `tests/` | gtest suites (`taupol_unit`) |

## Build

```bash
micromamba activate Nano
cd /data6/Users/youngwan/SKNanoAnalyzerV15
export SKNANO_ANALYSIS_MODULE_DIRS="/data6/Users/youngwan/Type1_Run3;/data6/Users/youngwan/WTauPolarization"
source setup.sh
./scripts/build.sh
```

List every module you still want installed: the manifest is rebuilt from the
modules in this variable.

## Test

```bash
ctest --test-dir "$SKNANO_BUILDDIR" -R taupol --output-on-failure
```

## Run

Tasks are selected with `--userflags`:

| flag | phase | status |
| --- | --- | --- |
| `GenStudy` (default) | generator level, MC only | implemented |

```bash
SKNano.py -a TauPolarization -i WtoLNu_2J -e 2024 -n 10 --reduction 100 --userflags GenStudy
```