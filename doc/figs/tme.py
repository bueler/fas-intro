#!/usr/bin/env python3

# generate a run-time figure based on results from fas/py/study/tme.sh

import sys
import numpy as np
import matplotlib.pyplot as plt

INFILE = 'results-tme.txt'

# columns: m derr err1 err2 err3
# where derr = discretization error |u-u_ex|_2 (i.e. error from F-cycle + 7 V(1,1))
#       err1 = error |u-u_ex|_2 from  F(1,1)          <-- 9 WU
#       err2 = error |u-u_ex|_2 from  F(1,0)          <-- 5 WU
#       err3 = error |u-u_ex|_2 from  F(1,0) + Rin    <-- 5 WU

LABELS = ['F(1,1) = 9 WU', 'F(1,0) = 5 WU', 'F(1,0)+$R_{inj}$ = 5 WU']
SYMBOLS = ['s','o','o']
MFC = ['w','k','w']
MSIZE = [10.0,10.0,10.0]

SHOW = False

def writeout(outname):
    if SHOW:
        plt.show()
    else:
        print('writing file ' + outname)
        plt.savefig(outname,bbox_inches='tight')

print('reading ' + INFILE + ' ...')
v = np.array(np.loadtxt(INFILE))
#print(v)
m = v[:,0]
de = v[:,1]

plt.figure(figsize=(7,6))
for k in range(3):
    plt.semilogx(m,v[:,k+2] / de,'k' + SYMBOLS[k],
                 mfc=MFC[k],label=LABELS[k],ms=MSIZE[k])
ax = plt.gca()
ax.set_xlim([1.0e2,1.0e6])
ax.set_ylim([1.0,2.0])
plt.xlabel('$m$',fontsize=16.0)
plt.ylabel('error relative to discretization',fontsize=16.0)
plt.legend(loc='lower left',fontsize=14.0)
plt.minorticks_off()
plt.grid(True)
writeout('tme.pdf')
