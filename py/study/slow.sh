#!/bin/bash
set -e
set +x

# solver complexity (optimality) results for NGS only
# shows work units and time

# run to generate results shown in ../../doc/fas.pdf:
#   $ ./slow.sh &> results-slow.txt
# and append these to results-optimal.txt

# NGS only; number of cycles determined using trial and error so that numerical
#   error was within a factor of two of error from F-cycles (see optimal.sh)
/usr/bin/time -f "\t%e real" ../fas1.py -mms -ngsonly -K 3 -cycles 20
/usr/bin/time -f "\t%e real" ../fas1.py -mms -ngsonly -K 4 -cycles 100
/usr/bin/time -f "\t%e real" ../fas1.py -mms -ngsonly -K 5 -cycles 800
/usr/bin/time -f "\t%e real" ../fas1.py -mms -ngsonly -K 6 -cycles 6400
/usr/bin/time -f "\t%e real" ../fas1.py -mms -ngsonly -K 7 -cycles 51200

