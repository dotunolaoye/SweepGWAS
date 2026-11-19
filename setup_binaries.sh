#!/usr/bin/env bash
# Puts the Linux builds of PLINK 2 and GCTA used in the paper into bin/, so
# that the app can find them on a server where neither is installed. The
# binaries are not kept in the repository because of their size.
#
# PLINK 2.0 v2.0.0-a.7.7 (18 Sep 2026) is downloaded from its official
# archive. GCTA 1.94.4 cannot be downloaded by script from the GCTA website,
# so its Linux build is copied from GCTA_SOURCE; download
# gcta-1.94.4-linux-kernel-3-x86_64.zip from
# https://yanglab.westlake.edu.cn/software/gcta/ and point GCTA_SOURCE at the
# gcta64 file inside it.
#
# Usage: GCTA_SOURCE=/path/to/gcta64 bash setup_binaries.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin="$here/bin"
mkdir -p "$bin"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# the generic x86_64 build rather than AVX2, so it runs on any server CPU
curl -fsSL -o "$tmp/plink2.zip" https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260918.zip
unzip -o -q "$tmp/plink2.zip" plink2 intel-simplified-software-license.txt -d "$tmp"
cp "$tmp/plink2" "$bin/plink2"
cp "$tmp/intel-simplified-software-license.txt" "$bin/PLINK2_INTEL_LICENSE.txt"
curl -fsSL -o "$bin/PLINK2_COPYING.txt" https://raw.githubusercontent.com/chrchang/plink-ng/master/2.0/COPYING

: "${GCTA_SOURCE:?set GCTA_SOURCE to the Linux gcta64 binary (GCTA 1.94.4)}"
cp "$GCTA_SOURCE" "$bin/gcta64"
curl -fsSL -o "$bin/GCTA_LICENSE.txt" https://raw.githubusercontent.com/jianyangqt/gcta/master/LICENSE

chmod 755 "$bin/plink2" "$bin/gcta64"
"$bin/plink2" --version
"$bin/gcta64" 2>&1 | grep -m1 "version" || true
