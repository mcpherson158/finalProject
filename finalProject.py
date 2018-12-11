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
5]
https://stackoverflow.com/questions/13013781/how-to-draw-a-rectangle-over-a-spe
cific-region-in-a-matplotlib-graph
6]
https://stackoverflow.com/questions/6697259/interactive-matplotlib-plot-with-tw
o-sliders
7]
https://stackoverflow.com/questions/3881453/numpy-add-row-to-array
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as d
import time
from matplotlib.patches import Rectangle
from matplotlib.widgets import Slider

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
                    flag = 'black'
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


class rectangleObject:

    def __init__(self, dimentions):
        xDim, yDim = dimentions
        self.xDim = xDim
        self.yDim = yDim
        self.yDimLower = yDim


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
        xDim = rectOb.xDim
        yDim = rectOb.yDim
        yDimLower = rectOb.yDimLower

        xpos = self.x[0]
        ypos = self.x[1]

        r = self.diameter / 2
        if xpos + r >= xDim:
            if self.v[0] > 0:
                self.v[0] = -self.v[0]
        elif xpos - r <= -xDim:
            if self.v[0] < 0:
                self.v[0] = -self.v[0]
        if ypos + r >= yDim:
            if self.v[1] > 0:
                self.v[1] = -self.v[1]
        elif ypos - r <= -yDimLower:
            if self.v[1] < 0:
                self.v[1] = -self.v[1]

    def speed(self, v):
        self.v = v

    def positionIncrease(self, dx):
        self.x = self.x + dx

    def velocityIncrease(self, dv):
        self.v = self.v + dv


class trackedParticle(robovac):

    def __init__(self, x, v, step, type, steps):
        robovac.__init__(self, v, x, step, type)
        self.xRecord = np.empty((steps, 2), dtype=float)
        self.index = 0
        self.diameter = 1.6

    def move(self, i, p, pX, delt):
        self.x += self.v*delt
        self.xRecord[self.index] = self.x
        self.index += 1


def stumble(mode='double', reaction=True, seed=0, length=20, dimentions=(10, 10), piston=True, steps=500, track=True, baseline=None, displacement = 0.01):
    # global xDim
    # global yDim
    global rectOb
    rectOb = rectangleObject(dimentions)

    xDim, yDim = dimentions
    np.random.seed(0)
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

    if track:
        x = np.array((float(np.random.rand(1)+displacement),
                      float(np.random.rand(1)+displacement) ))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))

        tracked = trackedParticle(x, v, step, 'tracked', steps)
        p = np.append(p, tracked)
        pX = np.vstack([pX, tracked.x])
        pT = np.append(pT, tracked.type)
        length += 1

    fig = plt.figure()
    ax = plt.axes()

    point1, = ax.plot([], [], 'bo', markersize=8)
    point2, = ax.plot([], [], 'ro', markersize=8)
    point3, = ax.plot([], [], 'go', markersize=8)
    point4, = ax.plot([], [], 'ko', markersize=8)
    point5, = ax.plot([], [], 'yo', markersize=14)
    linea,   = ax.plot([], [], 'orange')
    lineb,  = ax.plot([], [], 'c:')

    ax.set_xlim(-xDim-1, xDim+1)
    ax.set_ylim(-yDim-1, yDim+1)
    plt.title('Stumbling through the dark...')
    rectCentX, rectCentY = -(xDim), -(yDim)
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((rectCentX, rectCentY), 2*xDim, 2*yDim,
                          fill=False))

    if piston:
        x = (-xDim, xDim)
        y = (yDim, yDim)
        [line] = ax.plot(x, y)
        # Define an axes area and draw a slider in it
        amp_slider_ax = fig.add_axes([0.25, 0.02, 0.55, 0.03])
        amp_slider = Slider(amp_slider_ax, 'Height', -yDim, yDim,
                            valinit=rectOb.yDim)

        # Define an action for modifying the line when any slider value changes
        def sliders_on_changed(val):
            line.set_ydata(amp_slider.val)
            fig.canvas.draw_idle()
            rectOb.yDim = val
        amp_slider.on_changed(sliders_on_changed)

    dt = 0.01
    abcde = 0
    while abcde < steps:
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

                blackIndex = pT == 'black'
                xvalBlack = pX[blackIndex, 0]
                yvalBlack = pX[blackIndex, 1]
                point4.set_data(xvalBlack, yvalBlack)

        if track:
            xvalTrack = pX[-1, 0]
            yvalTrack = pX[-1, 1]
            point5.set_data(xvalTrack, yvalTrack)
            recIndex = p[-1].index
            xRec = p[-1].xRecord[:recIndex, 0]
            yRec = p[-1].xRecord[:recIndex, 1]
            linea.set_data(xRec, yRec)

            if baseline != None:
                xbase, ybase = baseline
                lineb.set_data(xbase, ybase)

        for i in range(len(p)):

            p[i].move(i, p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX, p, reaction)
        for j in range(length):
            pX[j] = p[j].x
        plt.pause(.0001)
        abcde += 1

def stumbleCompare(mode='double', reaction=True, seed=0, length=20, dimentions=(10, 10), piston=True, steps=500, track=True):
    # global xDim
    # global yDim
    global rectOb
    rectOb = rectangleObject(dimentions)

    xDim, yDim = dimentions
    np.random.seed(0)
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

    if track:
        x = np.array((float(np.random.rand(1)),
                      float(np.random.rand(1))))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))

        tracked = trackedParticle(x, v, step, 'tracked', steps)
        p = np.append(p, tracked)
        pX = np.vstack([pX, tracked.x])
        pT = np.append(pT, tracked.type)
        length += 1



    ######################################################################################################################

    rectCentX, rectCentY = -(xDim), -(yDim)



    dt = 0.01
    abcde = 0
    while abcde < steps:
        for types in range(len(p)):
            pT[types] = p[types].type

        blueIndex = pT == 'blue'

        xvalBlue = pX[blueIndex, 0]
        yvalBlue = pX[blueIndex, 1]


        if mode == 'double':
            redIndex = pT == 'red'
            xvalRed = pX[redIndex, 0]
            yvalRed = pX[redIndex, 1]


            if reaction:
                greenIndex = pT == 'green'
                xvalGreen = pX[greenIndex, 0]
                yvalGreen = pX[greenIndex, 1]


                blackIndex = pT == 'black'
                xvalBlack = pX[blackIndex, 0]
                yvalBlack = pX[blackIndex, 1]


        if track:
            xvalTrack = pX[-1, 0]
            yvalTrack = pX[-1, 1]

            recIndex = p[-1].index
            xRec = p[-1].xRecord[:recIndex, 0]
            yRec = p[-1].xRecord[:recIndex, 1]


        for i in range(len(p)):

            p[i].move(i, p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX, p, reaction)
        for j in range(length):
            pX[j] = p[j].x

        abcde += 1


    return p[-1]

def presentation(mode='double', reaction=True, seed=0, length=20, dimentions=(10, 10), piston=True, steps=500, track=True, displacement=0.01, compare=True):
    if compare:
        pX = stumbleCompare(mode, reaction, seed, length, dimentions, piston, steps, track)
        xplot = pX.xRecord[:, 0]
        yplot = pX.xRecord[:, 1]
        baseline = (xplot, yplot)
    else:
        baseline = None
        displacement = 0
    stumble(mode, reaction, seed, length, dimentions, piston, steps, track, baseline, displacement)


presentation()
