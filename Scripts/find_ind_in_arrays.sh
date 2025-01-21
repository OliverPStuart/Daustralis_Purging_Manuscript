#!/bin/bash

# Define three arrays
LHISI_ARR=(C01210 C01213 C01215 C01218 C01219 C01222 C01223 C01226 C01227 C01231 C01232 C10133)
LHIP_ARR=(C01211 C01216 C01217 C01220 C01221 C01224 C01225 C01228 C01230 C01233 C01234 C10223)
IND_ARR=(C01210 C01211 C01213 C01215 C01217 C01219 C01220 C01222 C01223 C01224 C01225 C01227 C01230 C01231 C01232 C01234 C10223 PAUL4 VAN2)

# Function to check if a value is in an array
contains_element() {
    local element match="$1"
    shift
    for element; do
        [[ "$element" == "$match" ]] && return 0
    done
    return 1
}

# Loop through each element in array3
for ind in "${IND_ARR[@]}"; do
    if contains_element "$ind" "${LHISI_ARR[@]}"; then
        echo "$ind is in LHISI"
    elif contains_element "$ind" "${LHIP_ARR[@]}"; then
        echo "$ind is in LHIP"
    else
        echo "$ind is in WILD"
        # Command(s) for when item is in neither array
    fi
done
