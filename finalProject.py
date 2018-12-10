# -*- coding: utf-8 -*-
"""

Sources:
1]
Stop clipping: https://seanny1986.wordpress.com/2017/10/01/simulation-of-elasti
c-collisions-in-python/
2]
cdist and using array: https://stackoverflow.com/questions/28144084/many-partic
les-in-box-physics-simulation
3]
Inelastic Collsions: http://www.batesville.k12.in.us/physics/APPhyNet/Dynamics/
Collisions/inelastic_2d.htm
4]
https://stackoverflow.com/questions/16094563/numpy-get-index-where-value-is-tru
e#16094877
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


def newParticle(p1, p2, v, flag):
    p1.v = v

    p1.type = flag

    p2.type = None
    p2.v[0] = 0
    p2.v[1] = 0
    p2.x[0] = (np.random.random() + 500)*1000
    p2.x[0] = (np.random.random() + 500)*1000


def collisionDetection(pX, p, reaction):

    dist = d.cdist(pX, pX)

    i = 0
    j = 0

    for i in range(len(p)):
        for j in range(i):

            collideDist = p[i].diameter/2 + p[j].diameter/2
            if dist[i, j] < collideDist:

                if not reaction:
                    collisionCalculaiton(p[i], p[j], dist[i, j], collideDist)
                elif (p[i].type == 'blue' and p[j].type == 'red') or        \
                     (p[i].type == 'red' and p[j].type == 'blue'):
                    flag = 'green'
                    collisionReaction(p[i], p[j], i, j, p, pX, flag)
                elif (p[i].type == 'blue' and p[j].type == 'green') or      \
                     (p[i].type == 'green' and p[j].type == 'blue'):
                    flag = 'purple'
                    collisionReaction(p[i], p[j], i, j, p, pX, flag)
                else:
                    collisionCalculaiton(p[i], p[j], dist[i, j], collideDist)


def collisionReaction(p1, p2, i, j, p, pX, flag):

    m1 = p1.m
    m2 = p2.m

    v1 = p1.v
    v2 = p2.v

    vCMx = (m1*v1[0]+m2*v2[0])/(m1+m2)
    vCMy = (m1*v1[1]+m2*v2[1])/(m1+m2)

    vCM = np.array((vCMx, vCMy))

    newParticle(p1, p2, vCM, flag)


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

    def __init__(self, v, x, step, type):
        self.x = np.array(x)
        self.v = np.array(v)
        self.diameter = 1
        self.m = 1
        self.bumpFlag = False
        self.time = time.time()
        self.type = type

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


def stumble(mode='single', reaction=False, seed=0, length=6):
    np.random.seed(0)
    # lengthBlue = 3
    # lengthRed = 3
    if mode == 'double':
        lengthBlue = int(length/2)
        lengthRed = int(length-lengthBlue)
        length = lengthBlue + lengthRed
    else:
        lengthBlue = length-1

    step = 40
    # v = (5, 5)
    p = [0]*length
    pX = np.zeros((length, 2))
    pT = np.empty((length), dtype=object)

    for j in range(lengthBlue+1):
        x = np.array((float(np.random.rand(1)),
                      float(np.random.rand(1))))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))
        # create an instance of the drunkard class
        p[j] = robovac(v, x, step, 'blue')
        pX[j, 0] = p[j].x[0]
        pX[j, 1] = p[j].x[1]

    if mode == 'double':
        for k in range(lengthRed):
            k += j
            x = np.array((float(np.random.rand(1))*8,
                          float(np.random.rand(1))*8))
            v = np.array((float(np.random.rand(1))*50,
                          float(np.random.rand(1))*50))
            p[k] = robovac(v, x, step, 'red')
            pX[k, 0] = p[k].x[0]
            pX[k, 1] = p[k].x[1]

    # fig = plt.figure()
    ax = plt.axes()

    point1, = ax.plot([], [], 'bo', markersize=5)
    point2, = ax.plot([], [], 'ro', markersize=5)
    point3, = ax.plot([], [], 'go', markersize=5)
    point4, = ax.plot([], [], 'ko', markersize=5)

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    plt.title('Stumbling through the dark...')

    dt = 0.01
    abcde = 0
    while abcde < 600:
        for types in range(len(p)):
            pT[types] = p[types].type

        blueIndex = pT == 'blue'
        xvalBlue = pX[blueIndex, 0]
        yvalBlue = pX[blueIndex, 1]
        point1.set_data(xvalBlue, yvalBlue)

        if mode == 'double':
            redIndex = pT == 'red'
            xvalRed = pX[redIndex, 0]
            yvalRed = pX[redIndex, 1]
            point2.set_data(xvalRed, yvalRed)

            if reaction:
                greenIndex = pT == 'green'
                xvalGreen = pX[greenIndex, 0]
                yvalGreen = pX[greenIndex, 1]
                point3.set_data(xvalGreen, yvalGreen)

                purpleIndex = pT == 'purple'
                xvalPurple = pX[purpleIndex, 0]
                yvalPurple = pX[purpleIndex, 1]
                point4.set_data(xvalPurple, yvalPurple)

        for i in range(len(p)):

            p[i].move(i, p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX, p, reaction)
        for j in range(length):
            pX[j] = p[j].x
        plt.pause(.0001)
        abcde += 1


stumble(mode='single', reaction=False, length=40)
