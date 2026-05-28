#!/usr/bin/env bash
set -euo pipefail

R_BIN="$(command -v R)"

if [[ -z "${R_BIN}" ]]; then
  echo "R not found in PATH; run this script inside the pixi environment." >&2
  exit 1
fi

ENV_PREFIX="${PIXI_ENV_PREFIX:-$(cd "$(dirname "${R_BIN}")/.." && pwd)}"
R_LIB="${ENV_PREFIX}/lib/R/library"
POST_LINK="${ENV_PREFIX}/bin/.bioconductor-genomeinfodbdata-post-link.sh"
STAGING="${ENV_PREFIX}/share/genomeinfodbdata-1.2.13"
TARBALL="${STAGING}/\"GenomeInfoDbData_1.2.13.tar.gz\""

if [[ -d "${R_LIB}/GenomeInfoDbData" ]]; then
  echo "GenomeInfoDbData already present at ${R_LIB}/GenomeInfoDbData"
  exit 0
fi

if [[ ! -f "${POST_LINK}" ]]; then
  echo "Post-link script not found: ${POST_LINK}" >&2
  exit 1
fi

PREFIX="${ENV_PREFIX}" PATH="${ENV_PREFIX}/bin:${PATH}" bash "${POST_LINK}" || true

if [[ ! -f "${TARBALL}" ]]; then
  echo "Expected tarball not found after post-link attempt: ${TARBALL}" >&2
  exit 1
fi

"${ENV_PREFIX}/bin/R" CMD INSTALL --library="${R_LIB}" "${TARBALL}"

echo "GenomeInfoDbData installed into ${R_LIB}"
