#!/bin/bash
set -e
set +x

# solver complexity (textbook multigrid efficiency) results for F-cycles

# run to generate results shown in ../../doc/fas.pdf:
#   $ ./tme.sh &> results-tme.txt

# fastest to discretization error: one F-cycle
for LEV in 7 8 9 10 11 12 13 14 15 16 17 18; do
    # first "-cyclemax 8" run is only to capture discretization error
    for OPT in "-cyclemax 8" "-cyclemax 1" "-cyclemax 1 -up 0" "-cyclemax 1 -up 0 -R inj"; do
        /usr/bin/time -f "\t%e real" ../fas1.py -mms -fcycle -rtol 0 -K $LEV $OPT
    done
done
