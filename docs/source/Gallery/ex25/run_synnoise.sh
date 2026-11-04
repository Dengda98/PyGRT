#!/bin/bash

set -euo pipefail

# BEGIN SYNNOISE
grt synnoise -Mmodel4 -Csources.txt -Ustations.txt -Oambient_run -N0/10
# END SYNNOISE
