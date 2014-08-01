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
data = np.genfromtxt(fname='en.dat') 
#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.subplot(211)
pl.xlim(0,4)
pl.ylabel('Energy [hartree]')
pl.plot(data[:,0],data[:,2],'b--',linewidth=2,label='Potential')
pl.plot(data[:,0],data[:,3],'g-.',linewidth=2,label='Quantum Potential')
pl.plot(data[:,0],data[:,4],'k-',linewidth=2,label='Energy')
lg = pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)
lg.draw_frame(False)

pl.subplot(212)
#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
pl.xlim(0,4)
pl.xlabel('time [a.u.]')
pl.ylabel('Energy [hartree]')
pl.plot(data[:,0],data[:,1],'r--',linewidth=2,label='Kinetic')
pl.yscale('log')
lg = pl.legend()
lg.draw_frame(False)
#pl.title()

pl.savefig('energy.pdf')
pl.show() 

