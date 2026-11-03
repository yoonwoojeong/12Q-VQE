#!/bin/sh
# Fetch libcint (headers + shared library) from conda-forge into third_party/libcint.
# Ubuntu/Debian do not package libcint; this uses the conda-forge build 6.1.3
# (requires libgcc >= 14 at runtime, which the build tells us about).
set -eu
PKG=libcint-6.1.3-h64ba90f_2.conda
URL=https://conda.anaconda.org/conda-forge/linux-64/$PKG
DEST=third_party/libcint
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$DEST"
curl -sSL -o "$TMP/$PKG" "$URL"
(cd "$TMP" && unzip -oq "$PKG" && for f in pkg-*.tar.zst; do zstd -dc "$f" | tar -x; done)
rm -rf "$DEST/include" "$DEST/lib"
cp -r "$TMP/include" "$TMP/lib" "$DEST/"
echo "libcint installed in $DEST"
