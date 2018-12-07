# -*- coding: utf-8 -*-
"""

Sources:
Stop clipping: https://seanny1986.wordpress.com/2017/10/01/simulation-of-elasti
c-collisions-in-python/
cdist and using array: https://stackoverflow.com/questions/28144084/many-partic
les-in-box-physics-simulation
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


def collisionDetection(pX, p):

    dist = d.cdist(pX, pX)

    for i in range(len(pX)):
        for j in range(i):
            collideDist = p[i].diameter/2 + p[j].diameter/2
            if dist[i, j] < collideDist:
                # print(dist[i,j],i,j)
                collisionCalculaiton(p[i], p[j], dist[i, j], collideDist)


def collisionCalculaiton(p1, p2, dist, collideDist):
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
        delv1 = -2*p2.m/total_m*np.inner(p1.v-p2.v, p1.x-p2.x) /   \
            np.sum((p1.x-p2.x)**2)*(p1.x-p2.x)
        delv2 = -2*p1.m/total_m*np.inner(p2.v-p1.v, p2.x-p1.x) /    \
            np.sum((p2.x-p1.x)**2)*(p2.x-p1.x)
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

    def move(self, i, p, pX, delt):
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

    def speed(self, v):
        self.v = v

    def positionIncrease(self, dx):
        self.x = self.x + dx

    def velocityIncrease(self, dv):
        self.v = self.v + dv
# Main program


def stumble():
    length = 10
    step = 40
    # v = (5, 5)
    p = [0]*length
    pX = np.zeros((length, 2))
    for j in range(length):
        x = np.array((float(np.random.rand(1)), float(np.random.rand(1))))
        v = np.array((float(np.random.rand(1))*50, float(np.random.rand(1))*50))
        p[j] = robovac(v, x, step)  # create an instance of the drunkard class
        pX[j, 0] = p[j].x[0]
        pX[j, 1] = p[j].x[1]

    fig = plt.figure()  # can use facecolor = 'white' argument to change figure background
    ax = plt.axes()

    point1, = ax.plot([], [], 'bo', markersize=5)
    # point2, = ax.plot([], [], 'bo' , markersize=12)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    plt.title('Stumbling through the dark...')

    dt = 0.01
    abcde = 0
    while abcde < 600:
        xval = pX[:, 0]
        yval = pX[:, 1]
        point1.set_data(xval, yval)

        for i in range(len(p)):

            p[i].move(i, p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX, p)
        for j in range(length):
            pX[j] = p[j].x
        plt.pause(.0001)  # this causes a stupid matplotlib warning - ignore!
        abcde += 1


stumble()
