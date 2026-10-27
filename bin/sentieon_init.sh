#!/bin/bash

# Written by Ramprasad Neethiraj and released under the MIT license.
# See git repository (https://github.com/nf-core/raredisease) for full license text.

# Sentieon initialization script
# This script takes as input the name of a environment
# variable holding the Sentieon license encoded as Base64 text
set -eu
LICENSE_ENCODED="${!1}"
if [ "${#LICENSE_ENCODED}" -lt "1500" ]; then  # Sentieon License server
    export SENTIEON_LICENSE=$(echo -e "$LICENSE_ENCODED" | base64 -d)
else  # Localhost license file
    export SENTIEON_LICENSE=$(mktemp)
    echo -e "$LICENSE_ENCODED" | base64 -d > $SENTIEON_LICENSE
fi
