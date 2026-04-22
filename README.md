# NGS De Novo Assembly Pipeline

Teaching-oriented bacterial assembly and polishing workflows for Ubuntu/WSL.

## Repository Layout

```text
NGS---de-novo-pipeline/
├── config.sh
├── setup.sh
├── envs/
│   ├── polishing.yml
│   └── busco.yml
├── scripts/
│   ├── 00_download_data.sh
│   ├── 01_pipeline_hybrid.sh
│   └── 02_pipeline_polish.sh
├── data/
│   ├── raw/
│   ├── reference/
│   └── db/
└── results/
```

## Quick Start

```bash
git clone https://github.com/yerpopek-cmyk/NGS---de-novo-pipeline.git
cd NGS---de-novo-pipeline
bash setup.sh
```

Download reads before running the pipelines:

```bash
bash scripts/00_download_data.sh SRRxxxxxxxx
```

`SRR25745292` is only a teaching example. Students should normally replace it with their own accession.

After downloading, review `config.sh` and update any organism-specific settings you need:

- `NANOPORE` if you also have long reads
- `GENOME_SIZE`
- `BUSCO_LINEAGE`
- `GTF_FILE`
- `REGION_CONTIG`, `REGION_START`, `REGION_END` for the optional reporting step

Run the pipelines from the repository root:

```bash
bash scripts/01_pipeline_hybrid.sh
bash scripts/02_pipeline_polish.sh
```

## WSL Notes

- Use the Ubuntu/WSL terminal to run the pipeline.
- Keep the repository inside the Linux filesystem when possible, for example under `/home/<user>/`.
- Shell and YAML files are pinned to LF line endings through `.gitattributes` to prevent `\r` issues.

## What `setup.sh` does

- Normalizes line endings on shell and YAML files.
- Creates or updates `polishing_env`.
- Creates or updates `busco_env`.
- Downloads the configured BUSCO lineage into `data/db/busco_downloads`.

## Environment Model

- `polishing_env` contains assembly, alignment, polishing, and visualization tools used by both pipelines.
- `busco_env` contains BUSCO and its heavier gene-prediction dependencies, isolated from the rest of the workflow to avoid solver conflicts.

## Resume / Checkpoints

Each pipeline writes checkpoint files under `results/.checkpoints/`.

- If a step already produced its expected output, the pipeline skips it on the next run.
- To force a complete re-run:

```bash
FORCE_RERUN=true bash scripts/01_pipeline_hybrid.sh
FORCE_RERUN=true bash scripts/02_pipeline_polish.sh
```
