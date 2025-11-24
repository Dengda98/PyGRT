#!/bin/bash
# Run all tests

set -euo pipefail

dirs=$(find . -maxdepth 1 -mindepth 1 -type d | sort)

for dir in $dirs; do
    cd $dir > /dev/null
    for fname in $(ls *.sh); do
        echo "@ --------------- $dir, $fname -------------"
        bash $fname
    done
    cd - > /dev/null
done