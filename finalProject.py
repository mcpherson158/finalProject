# -*- coding: utf-8 -*-
"""

Created on Thu Nov 10 19:26:26 2016

@author: jrathman
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.widgets import Slider, Button, RadioButtons

# import matplotlib.animation as animation

def partDistance(particle1,particle2):
    dist = np.sqrt( (particle1.xpos - particle2.xpos)**2 + (particle1.ypos - particle2.ypos))
    print(type(dist))
    return dist



class robovac:

    def __init__(self, batteryLife,xpos):
        self.batteryLife = batteryLife*60
        self.xpos = xpos
        self.ypos = 0
        self.theta = np.pi / 4
        self.v = .5
        self.diameter = 1
        self.mass = 1

    def collisionDetection(self, i, p):

        for j in range(len(p)):
            if j != i:

                dist = partDistance(p[i], p[j])
                
                if dist < self.diameter:
                    r1 = np.array((self.xpos, self.ypos))
                    r2 = np.array((p[j].xpos, p[j].ypos))

                    m1 = self.mass
                    m2 = p[j].mass

                    v1 = self.v
                    v2 = p[j].v

                    r_rel = r1 -r2
                    v_rel = v1 - v2

                    v_cm = (m1*v1 + m2*v2) / (m1 + m2)

                    rr_rel = np.dot(r_rel,r_rel)
                    vr_rel = np.dot(v_rel,r_rel)
                    v_rel = 2 * r_rel * vr_rel / rr_rel - v_rel

                    self.v = v_cm + v_rel * m2 / (m1 + m2)
    def move(self , i,p):

        # #Velocity in feet per second
        #
        # if self.xpos >= 10:
        #     if self.ypos >= 10:
        #         thetaRange = (np.pi, 3*np.pi/2)
        #
        #     elif self.ypos <= -10:
        #         thetaRange = (np.pi/2 , np.pi)
        #     else:
        #         thetaRange = (np.pi/2 , 3*np.pi/2)
        # elif self.ypos >= 10:
        #     thetaRange = (np.pi , 2*np.pi)
        #
        # elif self.xpos <= -10:
        #     if self.ypos >= 10:
        #         thetaRange = (3*np.pi/2, 2*np.pi)
        #     elif self.ypos <= -10:
        #         thetaRange = (0,np.pi/2)
        #     else:
        #         thetaRange = (3*np.pi/2 , 5*np.pi/2)
        # elif self.ypos <= -10:
        #     thetaRange = (0 , np.pi)
        #
        #
        if abs(self.xpos) >= 10 or abs(self.ypos) >= 10:
            thetaHigh , thetaLow = thetaRange
            #the following comes from plugging in generic high and low values
            #into point-slop formula for a line
            self.theta = (thetaHigh - thetaLow) * np.random.random() + thetaLow
        self.collisionDetection(i,p)
        self.xpos = self.v*np.cos(self.theta) + self.xpos
        self.ypos = self.v*np.sin(self.theta) + self.ypos




    def speed(self, v):
        self.v = v
#Main program
def stumble(batteryLife = 300):
    length = 5
    v=.5
    p = [0]*length
    px = [0]*length
    py = [0]*length

    for j in range(length):
        p[j] = robovac(batteryLife,np.random.randint(0,10)) #create an instance of the drunkard class
        px[j] = p[j].xpos
        py[j] = p[j].ypos

    fig = plt.figure() #can use facecolor = 'white' argument to change figure background
    ax = plt.axes()
    #ax = plt.axes(xlim=(-80, 80), ylim = (-80, 80))
    point1, = ax.plot([], [], 'bo' , markersize=12)
    # point2, = ax.plot([], [], 'bo' , markersize=12)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    plt.title('Stumbling through the dark...')

    # px = [particle1.xpos,particle2.xpos]
    # py = [particle1.ypos,particle2.ypos]
    # p = (particle1,particle2)

    # Define an axes area and draw a slider in it
    amp_slider_ax = fig.add_axes([0.25, 0.15, 0.45, 0.03])
    amp_slider = Slider(amp_slider_ax, 'Amp', 0.1, 1, valinit=v)

    def sliders_on_changed(val):
        for i in range(len(p)):
            p[i].speed(amp_slider.val)


        #fig.canvas.draw_idle()
    amp_slider.on_changed(sliders_on_changed)


    while True:

        point1.set_data(px, py)
        # point2.set_data(particle2.xpos, particle2.ypos)

        for i in range(len(p)):

            p[i].move(i,p)
        for j in range(length):
            px[j] = p[j].xpos
            py[j] = p[j].ypos
        plt.pause( 0.001)  #this causes a stupid matplotlib warning - ignore!

stumble()
