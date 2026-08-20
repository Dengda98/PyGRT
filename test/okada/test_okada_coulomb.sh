#!/bin/bash

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REFERENCE_ROOT="$HERE/coulomb_reference"

labels=(kode_100 kode_200 kode_300 kode_400 kode_500)
inputs=(
    "Example-1.inp"
    "Example-6(200Kode).inp"
    "Example-7(200-300kode).inp"
    "Example-Kode400.inp"
    "Example-9(Kode500).inp"
)

for input in "${inputs[@]}"; do
    if [[ ! -f "$REFERENCE_ROOT/input/$input" ]]; then
        echo "Missing Coulomb input: $REFERENCE_ROOT/input/$input" >&2
        exit 1
    fi
done

WORK_DIR=$(mktemp -d "$HERE/.coulomb_okada.XXXXXX")
cleanup() {
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

extract_faults() {
    local input=$1
    local output=$2

    {
        printf '# extracted Coulomb finite-fault table\n'
        printf '# x-start y-start x-fin y-fin kode value1 value2 dip top bot\n'
        awk '
            /X-start/ {
                in_table = 1
                skip_placeholder = 1
                next
            }
            in_table && /Grid Parameters/ { exit }
            in_table && skip_placeholder {
                skip_placeholder = 0
                next
            }
            in_table && NF >= 11 && $1 ~ /^[0-9]+$/ { print }
        ' "$input"
    } > "$output"

    if [[ "$(wc -l < "$output")" -le 2 ]]; then
        echo "No finite-fault rows extracted from $input" >&2
        exit 1
    fi
}

extract_receivers() {
    local reference=$1
    local output=$2

    awk 'NR > 3 && NF >= 6 {
        printf "%.17g %.17g %.17g\n", $2, $1, -$3
    }' "$reference" > "$output"

    if [[ ! -s "$output" ]]; then
        echo "No receiver points extracted from $reference" >&2
        exit 1
    fi
}

for index in "${!labels[@]}"; do
    label=${labels[$index]}
    input=${inputs[$index]}
    reference_dir="$REFERENCE_ROOT/result/$label"
    fault_path="$WORK_DIR/faults_$label.inp"
    receiver_path="$WORK_DIR/receivers_$label.txt"
    output_dir="$WORK_DIR/$label"
    mkdir -p "$output_dir"

    extract_faults "$REFERENCE_ROOT/input/$input" "$fault_path"
    extract_receivers "$reference_dir/Displacement.cou" "$receiver_path"

    grt okada -I6/3.4641016151377544/2.7 -C"$fault_path" \
        -Q"$receiver_path" -N -e -O"$output_dir/pygrt.nc" -s
done

python -u "$HERE/test_okada_coulomb.py" "$REFERENCE_ROOT" "$WORK_DIR"
echo "test_okada_coulomb.sh: all checks passed"
