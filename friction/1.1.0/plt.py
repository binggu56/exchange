##!/usr/bin/python

import numpy as np
import pylab as pl
import seaborn as sns

sns.set_context("poster")
sns.set_style("ticks") # darkgrid, whitegrid, dark, white, and ticks
#with open("traj.dat") as f:
#    data = f.read()
#
#    data = data.split('\n')
#
#    x = [row.split(' ')[0] for row in data]
#    y = [row.split(' ')[1] for row in data]
#
#    fig = plt.figure()
#
#    ax1 = fig.add_subplot(111)
#
#    ax1.set_title("Plot title...")    
#    ax1.set_xlabel('your x label..')
#    ax1.set_ylabel('your y label...')
#
#    ax1.plot(x,y, c='r', label='the data')
#
#    leg = ax1.legend()
#fig = plt.figure()
data = np.genfromtxt(fname='traj.dat') 
#data = np.loadtxt('traj.dat')
for x in range(1,20):
    pl.plot(data[:,0],data[:,x])

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.legend()
pl.xlabel('time [a.u.]')
pl.ylabel('Quantum trajectories in x axis [bohr]')
pl.savefig('traj_lqf.pdf')
pl.show() 

