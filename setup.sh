#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"
CONFIG="${PROJECT_ROOT}/config.sh"

die() {
    echo "ERROR: $1" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

if [ -f "$CONFIG" ]; then
    # shellcheck source=/dev/null
    source "$CONFIG"
fi

BUSCO_ENV="${BUSCO_ENV:-busco_env}"
POLISH_ENV="${POLISH_ENV:-polishing_env}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-vibrio_odb12}"
BUSCO_DB_DIR="${BUSCO_DB_DIR:-data/db/busco_downloads}"

case "$BUSCO_DB_DIR" in
    /*) ;;
    *) BUSCO_DB_DIR="${PROJECT_ROOT}/${BUSCO_DB_DIR}" ;;
esac

require_command micromamba
require_command git

mkdir -p "${PROJECT_ROOT}/data/db" "${PROJECT_ROOT}/results"

normalize_lf() {
    local file="$1"
    if [ -f "$file" ]; then
        sed -i 's/\r$//' "$file"
    fi
}

echo "Normalizing line endings for shell and YAML files..."
while IFS= read -r file; do
    normalize_lf "$file"
done < <(find "$PROJECT_ROOT" -maxdepth 2 -type f \( -name "*.sh" -o -name "*.yml" -o -name "*.yaml" \))

create_or_update_env() {
    local env_name="$1"
    local env_file="$2"

    if micromamba env list | awk '{print $1}' | grep -Fxq "$env_name"; then
        echo "Updating existing environment: $env_name"
        micromamba install -y -n "$env_name" -f "$env_file"
    else
        echo "Creating environment: $env_name"
        micromamba create -y -n "$env_name" -f "$env_file"
    fi
}

create_or_update_env "$POLISH_ENV" "${PROJECT_ROOT}/envs/polishing.yml"
create_or_update_env "$BUSCO_ENV" "${PROJECT_ROOT}/envs/busco.yml"

mkdir -p "$BUSCO_DB_DIR"

if [ ! -d "${BUSCO_DB_DIR}/lineages/${BUSCO_LINEAGE}" ]; then
    echo "Downloading BUSCO lineage: ${BUSCO_LINEAGE}"
    micromamba run -n "$BUSCO_ENV" busco \
        --download "$BUSCO_LINEAGE" \
        --download_path "$BUSCO_DB_DIR"
else
    echo "BUSCO lineage already present: ${BUSCO_DB_DIR}/lineages/${BUSCO_LINEAGE}"
fi

echo ""
echo "Setup complete."
echo "  PROJECT_ROOT = ${PROJECT_ROOT}"
echo "  POLISH_ENV   = ${POLISH_ENV}"
echo "  BUSCO_ENV    = ${BUSCO_ENV}"
echo "  BUSCO_DB_DIR = ${BUSCO_DB_DIR}"
echo ""
echo "Next steps:"
echo "  bash scripts/01_pipeline_hybrid.sh"
echo "  bash scripts/02_pipeline_polish.sh"
