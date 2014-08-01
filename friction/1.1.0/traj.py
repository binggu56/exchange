##!/usr/bin/python

import numpy as np
import pylab as pl

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
pl.subplot(211)
for x in range(1,6):
    pl.plot(data[:,0],data[:,x],'k--')

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.ylim(0,1)
pl.xlabel('time [a.u.]')
pl.ylabel('Positions for $x$ [bohr]')
#pl.title('traj')
pl.subplot(212)
for x in range(7,12):
    pl.plot(data[:,0],data[:,x],'k--')

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.ylim(0,1)
pl.xlabel('time [a.u.]')
pl.ylabel('Positions for $x$ [bohr]')
pl.savefig('traj.pdf')
pl.show() 

