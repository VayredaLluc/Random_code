#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 14:42:17 2024

@author: lluc
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

data = np.loadtxt('double_pendulum_res_1.dat')
Nt = np.shape(data)[0]

print(Nt)
x1 = np.zeros((10,Nt))
y1 = np.zeros((10,Nt))
x2 = np.zeros((10,Nt))
y2 = np.zeros((10,Nt))


for i in range(10):
    data = np.loadtxt('double_pendulum_res_'+str(i+1)+'.dat')
    x1[i,:] = data[:,1]
    y1[i,:] = data[:,2]
    x2[i,:] = data[:,3]
    y2[i,:] = data[:,4]
    t = data[:,0]

def update(frame):
    k=frame*2
    plt.cla()
    
    ini = 0
    trans = 0.5
    if k <= 30:
        plt.plot(x1[ini,k],y1[ini,k],'r-',alpha=trans)
        plt.plot(x2[ini,0:k],y2[ini,0:k],'r-',alpha=trans)
    else:
        plt.plot(x1[ini,k],y1[ini,k],'r-',alpha=trans)
        plt.plot(x2[ini,(k-15):k],y2[ini,(k-15):k],'r-',alpha=trans)
        plt.plot(x2[ini,0:k],y2[ini,0:k],'r-',alpha=trans**4)
    
    x_line_12 = [x1[ini,k], x2[ini,k]]
    y_line_12 = [y1[ini,k], y2[ini,k]]

    x_line_01 = [x1[ini,k], 0]
    y_line_02 = [y1[ini,k], 0]

    plt.plot(x_line_12, y_line_12,'r-',alpha=trans)
    plt.plot(x_line_01, y_line_02,'r-',alpha=trans)
    plt.plot(x1[ini,k],y1[ini,k],'ro',alpha=trans)
    plt.plot(x2[ini,k],y2[ini,k],'ro',alpha=trans)
    plt.axis([-2.2,2.2,-2.2,2.2])
    
    ini = 1
    if k <= 15:
        plt.plot(x1[ini,k],y1[ini,k],'b-',alpha=trans)
        plt.plot(x2[ini,0:k],y2[ini,0:k],'b-',alpha=trans)
    else:
        plt.plot(x1[ini,k],y1[ini,k],'b-',alpha=trans)
        plt.plot(x2[ini,(k-15):k],y2[ini,(k-15):k],'b-',alpha=trans)
        plt.plot(x2[ini,0:k],y2[ini,0:k],'b-',alpha=trans**4)
    
    x_line_12 = [x1[ini,k], x2[ini,k]]
    y_line_12 = [y1[ini,k], y2[ini,k]]

    x_line_01 = [x1[ini,k], 0]
    y_line_02 = [y1[ini,k], 0]

    plt.plot(x_line_12, y_line_12,'b-',alpha=trans)
    plt.plot(x_line_01, y_line_02,'b-',alpha=trans)
    plt.plot(x1[ini,k],y1[ini,k],'bo',alpha=trans)
    plt.plot(x2[ini,k],y2[ini,k],'bo',alpha=trans)
    plt.axis([-2.2,2.2,-2.2,2.2])

    

fig = plt.figure(figsize=(8, 8))
ax1 = plt.subplot()

Writer = ani.writers['ffmpeg']
writer = Writer(fps=32, metadata=dict(artist='Me'), bitrate=1800)

Nframes = 10000
anim = ani.FuncAnimation(fig, update, 
                               frames = int((Nframes-1)/2), 
                               blit = False, interval=50)

anim.save('double_pendulum_close_dif.mp4', writer=writer)


    
    

