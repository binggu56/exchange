#!/usr/bin/python

import numpy as np
import pylab as pl
#import seaborn as sns

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
data = np.genfromtxt(fname='en.dat')
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.ylabel('Energy [hartree]')
pl.plot(data[:,0],data[:,2],'b--',linewidth=2,label='Potential')
pl.plot(data[:,0],data[:,3],'g-.',linewidth=2,label='Quantum Potential')
pl.plot(data[:,0],data[:,4],'k-',linewidth=2,label='Energy')
pl.plot(data[:,0],data[:,1],'r--',linewidth=2,label='Kinetic')
pl.legend()

pl.savefig('energy.pdf')
pl.show() 

