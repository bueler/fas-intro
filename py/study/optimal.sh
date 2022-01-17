#!/bin/bash
set -e
set +x

# solver complexity (optimality) results for F-cycles and V-cycles
# results graph shows time-vs-m

# run to generate results shown in ../../doc/fas.pdf:
#   $ ./optimal.sh &> results-optimal.txt

# fastest to discretization error: F-cycles with 3 V(1,0) cycles
for LEV in 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18; do
    /usr/bin/time -f "\t%e real" ../fas1.py -mms -fcycle -cycles 3 -up 0 -K $LEV
done

# 12 V(1,1) cycles get discretization error too
for LEV in 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
    /usr/bin/time -f "\t%e real" ../fas1.py -mms -cycles 12 -K $LEV
done

