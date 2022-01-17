#!/usr/bin/env python3

# generate a run-time figure based on results from fas/py/study/{optimal.sh|slow.sh}

import numpy as np
import matplotlib.pyplot as plt
import sys

INFILE = 'results-optimal.txt'
# columns: method m WU l2err time
# method=0:  F-cycles with down=1,up=0 followed by 3xV(1,0)
# method=1:  12 V(1,1)-cycles
# method=2:  NGS with enough cycles to get 2*(discretization error)

LABELS = ['F-cycles+3xV(1,0)', '12 V(1,1) cycles', 'NGS sweeps']
SYMBOLS = ['s','.','o']
MFC = ['k','k','w']
MSIZE = [9.0,16.0,10.0]

SHOW = False

def writeout(outname):
    if SHOW:
        plt.show()
    else:
        print('writing file ' + outname)
        plt.savefig(outname,bbox_inches='tight')

print('reading ' + INFILE + ' ...')
v = np.reshape(np.array(np.loadtxt(INFILE)), (-1,5))
#print(v)

method = v[:,0]
m = v[:,1]
h = 1.0 / m
wu = v[:,2]
err = v[:,3]
time = v[:,4]

plt.figure(figsize=(7,6))
for k in range(3):
    plt.loglog(m[method==k],time[method==k],'k' + SYMBOLS[k],
               mfc=MFC[k],label=LABELS[k],ms=MSIZE[k])
    if k < 2:
        mm = m[method==k]
        ttime = time[method==k]
        tbig = ttime[ttime>2.0]
        mbig = mm[ttime>2.0]
        p = np.polyfit(np.log(mbig),np.log(tbig),1)
        plt.loglog(mbig,np.exp(p[0]*np.log(mbig)+p[1]),'k:')
        if k == 0:
            plt.text(8.0e4,4.0,r'$O(m^{%.2f})$' % p[0],fontsize=16.0)
        elif k == 1:
            plt.text(5.0e3,30.0,r'$O(m^{%.2f})$' % p[0],fontsize=16.0)

ax = plt.gca()
ax.set_xlim([10.0,1.0e6])
ax.set_ylim([0.1,140.0])
plt.xticks([10,100,1.0e3,1.0e4,1.0e5,1.0e6],fontsize=14.0)
plt.yticks([0.1,1.0,10.0,100.0],fontsize=12.0)
plt.gca().set_yticklabels(['0.1','1','10','100'],fontsize=14.0)
plt.xlabel('$m$',fontsize=16.0)
plt.ylabel('time (seconds)',fontsize=16.0)
plt.legend(loc='best',fontsize=14.0)
plt.minorticks_off()
plt.grid(True)

writeout('optimal.pdf')

