#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

mkdir -p "${ROOT_DIR}/data/raw/remote" "${ROOT_DIR}/data/raw/supplements"

fetch_file() {
  local url="$1"
  local out="$2"
  if [ -s "${out}" ]; then
    echo "[fetch_real_data] exists: ${out}"
    return
  fi
  curl --fail --location --retry 3 --output "${out}" "${url}"
  echo "[fetch_real_data] downloaded: ${out}"
}

verify_sha256() {
  local expected="$1"
  local file="$2"
  local observed
  observed="$(sha256sum "${file}" | awk '{print $1}')"
  if [ "${observed}" != "${expected}" ]; then
    echo "[fetch_real_data] checksum mismatch: ${file}" >&2
    echo "  expected: ${expected}" >&2
    echo "  observed: ${observed}" >&2
    return 1
  fi
  echo "[fetch_real_data] checksum verified: ${file}"
}

fetch_file \
  "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02811-x" \
  "${ROOT_DIR}/data/raw/remote/paper_article.html"
fetch_file \
  "https://api.github.com/repos/vagarwal87/saluki_paper" \
  "${ROOT_DIR}/data/raw/remote/github_repo.json"
fetch_file \
  "https://zenodo.org/api/records/7158836" \
  "${ROOT_DIR}/data/raw/remote/zenodo_record.json"
fetch_file \
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02811-x/MediaObjects/13059_2022_2811_MOESM2_ESM.xlsx" \
  "${ROOT_DIR}/data/raw/supplements/13059_2022_2811_MOESM2_ESM.xlsx"
fetch_file \
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-022-02811-x/MediaObjects/13059_2022_2811_MOESM3_ESM.xlsx" \
  "${ROOT_DIR}/data/raw/supplements/13059_2022_2811_MOESM3_ESM.xlsx"

verify_sha256 \
  "ddf9441d1fe85d091d4b7f0a444b8dcf3658aa1d4158f66948bff00d5406de98" \
  "${ROOT_DIR}/data/raw/supplements/13059_2022_2811_MOESM2_ESM.xlsx"
verify_sha256 \
  "00c49a9869e982527b07ac27896340fa34a3188199551f386027a8e8eb4f2217" \
  "${ROOT_DIR}/data/raw/supplements/13059_2022_2811_MOESM3_ESM.xlsx"
