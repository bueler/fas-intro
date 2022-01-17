#!/bin/bash
set -e
set +x

# convergence results in -mms case

# run to generate results shown in ../../doc/fas.pdf:
#   $ ./converge.sh &> converge.txt

# brutal NGS sweeps
for LEV in 3 4 5 6 7 8; do
    ../fas1.py -cycles 10000 -ngsonly -mms -K $LEV
done

# default V-cycles
for LEV in 3 4 5 6 7 8 9 10 11 12 13 14; do
    ../fas1.py -mms -cycles 12 -K $LEV
done

