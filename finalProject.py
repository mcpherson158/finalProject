# -*- coding: utf-8 -*-
"""

Sources:
Stop clipping: https://seanny1986.wordpress.com/2017/10/01/simulation-of-elastic-collisions-in-python/
cdist and using array: https://stackoverflow.com/questions/28144084/many-particles-in-box-physics-simulation
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as d
import time

from matplotlib.widgets import Slider, Button, RadioButtons

# import matplotlib.animation as animation


def partDistance(particle1, particle2):
    dist = np.sqrt((particle1.x[0] - particle2.x[0])**2 + (particle1.x[1] -
                                                           particle2.x[1])**2)

    return dist

def magnitude(vector):
    x, y = vector
    magnitude = np.sqrt(x**2 + y**2)
    return magnitude

def collisionDetection(pX,p):

    dist = d.cdist(pX, pX)

    for i in range(len(pX)):
        for j in range(i):
            collideDist = p[i].diameter/2 + p[j].diameter/2
            if dist[i,j] < collideDist:
                #print(dist[i,j],i,j)
                collisionCalculaiton(p[i],p[j], dist[i,j], collideDist)




def collisionCalculaiton(p1, p2, dist, collideDist):
    #
    # r1 = p1.diameter / 2
    # r2 = p2.diameter / 2
    # v1 = p1.v
    # v2 = p2.v
    # m1 = p1.m
    # m2 = p2.m
    #
    #
    # x1 = np.array(p1.x)
    # x2 = np.array(p2.x)
    # if np.allclose(x1,x2):
    #     x1 += (v1/magnitude(v1))*0.00001
    #
    #
    #
    # massTerm1 = (2*m2)/(m1+m2)
    # massTerm2 = (2*m1)/(m1+m2)
    # vsub1 = v1-v2
    # vsub2 = v2-v1
    # xsub1 = x1-x2
    # xsub2 = x2-x1
    # xmag1 = magnitude(xsub1)
    # xmag2 = magnitude(xsub2)
    # xmagsqr1 = xmag1**2
    # xmagsqr2 = xmag2**2
    # dotprod1 = np.dot(vsub1,xsub1)
    # dotprod2 = np.dot(vsub2,xsub2)
    #
    # prodTerm1 = dotprod1/xmagsqr1
    # prodTerm2 = dotprod2/xmagsqr2
    #
    #
    # v1 = v1 - massTerm1*prodTerm1*xsub1
    # v2 = v2 - massTerm2*prodTerm2*xsub2
    # print(v1)
    #
    # # if p1.bumpFlag == False:
    #
    #
    # p1.v = v1
    #     # p1.v =  v1 +v2
    #     # p1.bumpFlag = True
    #     # p1.time = time.time()
    # # if p2.bumpFlag == False:
    # p2.v = v2
    #     # p2.v = v2 - massTerm2*prodTerm2*xsub2
    #
    #     #p2.bumpFlag = True
    #     # p2.time = time.time()
    # # else:
    #     #p1.bumpFlag = False
    #     #p2.bumpFlag = False
    # xpos1 = p1.v[0] + p1.x[0]
    # ypos1 = p1.v[1] + p1.x[1]
    #
    # xpos2 = p2.v[0] + p2.x[0]
    # ypos2 = p2.v[1] + p2.x[1]
    #
    # difference = np.array(x1) - np.array(x2)
    # mag = magnitude(difference)
    # offset = mag - (p1.diameter/2 + p2.diameter/2)
    #
    # percentage1 = magnitude(p1.v) / (magnitude(p1.v) + magnitude(p2.v))
    #
    # p1.positionIncrease((-difference/mag)*offset * magnitude(p1.v)/2)
    # p2.positionIncrease((difference/mag)*offset * magnitude(p2.v)/2)
    # # if np.sqrt((xpos1 - xpos2)**2 + (ypos1 - ypos2)**2) < collideDist:
    #     #p1.x = p1.x + p1.v*collideDist/2
    #     #p2.x = p2.x + p2.v*collideDist/2
    #     #print(p1.v)

    delx = p1.x - p2.x
    distance = np.sqrt(np.sum(delx**2))
    if distance < p1.diameter/2 + p2.diameter/2:
        if distance == 0:
            p1.x -= 0.00001

        # orignally, when using the wikipedia equaiton, I missread the equaiton
        # and made the v-prime the new v, when it should have been deltaVself.
        # Source[ ]  showed me this mistake
        offset = dist - (p1.diameter/2 + p2.diameter/2)
        p1.positionIncrease((-delx/distance)*offset/2)
        p2.positionIncrease((delx/distance)*offset/2)
        total_m = p1.m + p2.m
        delv1 = -2*p2.m/total_m*np.inner(p1.v-p2.v,p1.x-p2.x)/np.sum((p1.x-p2.x)**2)*(p1.x-p2.x)
        delv2 = -2*p1.m/total_m*np.inner(p2.v-p1.v,p2.x-p1.x)/np.sum((p2.x-p1.x)**2)*(p2.x-p1.x)
        p1.velocityIncrease(delv1)
        p2.velocityIncrease(delv2)

class robovac:

    def __init__(self, v, x, step):
        self.x = np.array(x)
        self.v = np.array(v)
        self.diameter = 1
        self.m = 1
        self.bumpFlag = False
        self.time = time.time()




    def move(self, i, p,pX, delt):
        #xpos = self.x[0]
        #ypos = self.x[1]
        #xpos = self.v[0]*delt + xpos
        #ypos = self.v[1]*delt + ypos

        self.x += self.v*delt

    def wallCollision(self):
        xpos = self.x[0]
        ypos = self.x[1]

        r = self.diameter / 2
        if xpos + r >= 10:
            if self.v[0] > 0:
                self.v[0] = -self.v[0]
        elif xpos - r <= -10:
            if self.v[0] < 0:
                self.v[0] = -self.v[0]
        if ypos + r >= 10:
            if self.v[1] > 0:
                self.v[1] = -self.v[1]
        elif ypos - r <= -10:
            if self.v[1] < 0:
                self.v[1] = -self.v[1]


        # if abs(xpos) >= 10 or abs(ypos) >= 10:
        #     thetaHigh, thetaLow = thetaRange
        #     # the following comes from plugging in generic high and low values
        #     # into point-slop formula for a line
        #     theta = (thetaHigh - thetaLow) * np.random.random() + thetaLow
        # self.collisionDetection(i,p)
        #self.v[0] = speed*np.cos(theta)
        #self.v[1] = speed*np.sin(theta)
        # self.v = np.array(v0, v1)
        # self.x[0] = self.v[0] + xpos
        # self.x[1] = self.v[1] + ypos





        # collisionDetection(pX, p)

    def speed(self, v):
        self.v = v

    def positionIncrease(self, dx):
        self.x = self.x + dx

    def velocityIncrease(self,dv):
        self.v = self.v + dv
# Main program


def stumble():
    length = 100
    step=40
    # v = (5, 5)
    p = [0]*length
    pX = np.zeros((length,2))
    for j in range(length):
        x = np.array(( float(np.random.rand(1)),float(np.random.rand(1))))
        v = np.array(( float(np.random.rand(1))*50,float(np.random.rand(1))*50))
        p[j] = robovac(v, x, step) # create an instance of the drunkard class
        pX[j,0] = p[j].x[0]
        pX[j,1] = p[j].x[1]

    #p = [0]*2
    #pX = np.zeros((2,2))

    # p[0]=robovac((15,0),(0,0),step)
    # p[1]=robovac((0,0),(1,0),step)
    # for j in range(len(p)):
    #     pX[j,0] = p[j].x[0]
    #     pX[j,1] = p[j].x[1]
    #
    # print(pX)


    fig = plt.figure() #can use facecolor = 'white' argument to change figure background
    ax = plt.axes()
    #ax = plt.axes(xlim=(-80, 80), ylim = (-80, 80))
    point1, = ax.plot([], [], 'bo' , markersize=5)
    # point2, = ax.plot([], [], 'bo' , markersize=12)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    plt.title('Stumbling through the dark...')

    # px = [particle1.xpos,particle2.xpos]
    # py = [particle1.ypos,particle2.ypos]
    # p = (particle1,particle2)

    # Define an axes area and draw a slider in it
    # amp_slider_ax = fig.add_axes([0.25, 0.15, 0.45, 0.03])
    # amp_slider = Slider(amp_slider_ax, 'Amp', 0.1, 1, valinit=v)
    #
    # def sliders_on_changed(val):
    #     for i in range(len(p)):
    #         p[i].speed(amp_slider.val)
    #
    #
    #     #fig.canvas.draw_idle()
    # amp_slider.on_changed(sliders_on_changed)
    dt = 0.01
    abcde = 0
    while abcde < 600:
        xval = pX[:,0]
        yval = pX[:,1]
        point1.set_data(xval,yval)
        # point2.set_data(particle2.xpos, particle2.ypos)

        for i in range(len(p)):

            p[i].move(i,p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX,p)
        for j in range(length):
            pX[j] = p[j].x
        plt.pause(.0000001)  #this causes a stupid matplotlib warning - ignore!
        abcde += 1
stumble()
