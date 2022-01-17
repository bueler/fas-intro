#!/usr/bin/env python3

# generate a convergence figure based on results from fas/py/study/converge.sh

import numpy as np
import matplotlib.pyplot as plt
import sys

INFILE = 'results-converge.txt'
LABELS = [r'$10^4$ NGS sweeps', '12 V(1,1) cycles']
SYMBOLS = ['o','.']
MFC = ['w','k']
MSIZE = [14.0,14.0]

SHOW = False

def writeout(outname):
    if SHOW:
        plt.show()
    else:
        print('writing file ' + outname)
        plt.savefig(outname,bbox_inches='tight')

print('reading ' + INFILE + ' ...')
v = np.reshape(np.array(np.loadtxt(INFILE)), (-1,3))
# columns of v are:   method, m, |u-u_ex|_2
#print(v)

method = v[:,0]
m = v[:,1]
h = 1.0 / m
err = v[:,2]

plt.figure(figsize=(7,6))
for k in range(2):
    plt.loglog(h[method==k],err[method==k],'k' + SYMBOLS[k],
               mfc=MFC[k],label=LABELS[k],ms=MSIZE[k])
p = np.polyfit(np.log(h[method==1]),np.log(err[method==1]),1)
plt.loglog(h[method==1],np.exp(p[0]*np.log(h[method==1])+p[1]),'k:')
plt.text(4.0e-4,2.0e-7,'$O(h^{%.2f})$' % p[0],fontsize=16.0)

ax = plt.gca()
ax.set_xlim([1.0e-5,1.0e-1])
ax.set_ylim([1.0e-9,1.0e-1])
plt.xticks([1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1],fontsize=14.0)
plt.xlabel('$h$',fontsize=16.0)
plt.yticks([1.0e-9,1.0e-7,1.0e-5,1.0e-3,1.0e-1],fontsize=14.0)
plt.ylabel('numerical error',fontsize=16.0)
plt.minorticks_off()
plt.grid(True)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='best',fontsize=14.0)


writeout('converge.pdf')

