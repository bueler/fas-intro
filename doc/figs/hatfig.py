#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"

fsize=16.0
bigfsize=20.0

debug = False
def figsave(name):
    if debug:
        plt.show()  # debug
    else:
        print('saving %s ...' % name)
        plt.savefig(name,bbox_inches='tight',transparent=True)

def meshaxes(m):
    x = np.linspace(0.0,1.0,2)
    y = np.linspace(0.0,1.0,2)
    plt.plot(x, np.zeros(2), 'k', lw=1.0)
    plt.plot(np.zeros(2), y, 'k', lw=1.0)
    plt.text(-0.01,-0.2,'0',fontsize=fsize)
    plt.text(0.99,-0.2,'1',fontsize=fsize)
    xdots = np.linspace(0.0,1.0,m+1)
    return xdots

def hats(xdots,kk,style='k',width=2.0,dotstyle='k.',dotsize=12.0):
    for k in kk:
        plt.plot(xdots[k:k+3],[0.0,1.0,0.0],style,lw=width)
    plt.plot(xdots[1:-1],np.zeros(len(xdots)-2),dotstyle,ms=dotsize)
    plt.plot([0.0,1.0],[0.0,0.0],dotstyle,ms=dotsize,mfc='w')

m = 8
kk = [2, 3, 4]

# one hat figure
plt.figure(figsize=(10,4))
xdots = meshaxes(m)
hats(xdots,[2,])
plt.text(xdots[3]-0.01,-0.15,r'$x_p$',fontsize=bigfsize)
plt.text(xdots[3]+0.06,0.7,r'$\psi_p(x)$',fontsize=bigfsize)
plt.text(-0.05,0.95,'1',fontsize=fsize)
plt.axis([-0.2,1.2,-0.4,1.2])
plt.axis('off')
figsave('onehat.pdf')

# hat combination figure
plt.figure(figsize=(10,4))
xdots = meshaxes(m)
hats(xdots,kk,style='k:')
hats(xdots[::2],[1,])
plt.text(xdots[3]-0.06,-0.2,r'$2q-1$',fontsize=fsize)
plt.text(xdots[4]-0.02,-0.2,r'$2q$',fontsize=fsize)
plt.text(xdots[5]-0.05,-0.2,r'$2q+1$',fontsize=fsize)
plt.text(-0.05,0.95,'1',fontsize=fsize)
plt.axis([-0.2,1.2,-0.4,1.2])
plt.axis('off')
figsave('hatcombination.pdf')

