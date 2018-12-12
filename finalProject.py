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
Index for printing
https://stackoverflow.com/questions/16094563/numpy-get-index-where-value-is-tru
e#16094877
5]
Drawing a rectangle
https://stackoverflow.com/questions/13013781/how-to-draw-a-rectangle-over-a-spe
cific-region-in-a-matplotlib-graph
6]
Interactive sliders
https://stackoverflow.com/questions/6697259/interactive-matplotlib-plot-with-tw
o-sliders
7]
Add row to numpy array
https://stackoverflow.com/questions/3881453/numpy-add-row-to-array
8]
Math for elactic collisions
https://en.wikipedia.org/wiki/Elastic_collision
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as d
import time
from matplotlib.patches import Rectangle
from matplotlib.widgets import Slider


def partDistance(particle1, particle2):
    '''
    calculates distance between two particles
    '''
    dist = np.sqrt((particle1.x[0] - particle2.x[0])**2 + (particle1.x[1] -
                                                           particle2.x[1])**2)

    return dist


def magnitude(vector):
    '''
    calculates the magnitude of a vector
    '''
    x, y = vector
    magnitude = np.sqrt(x**2 + y**2)
    return magnitude


def newParticle(p1, p2, v, flag):
    '''
    sets one particle to the new one and the other to none type after a
    reaction occurs. Simulates two particles reacting to form a single one
    '''

    # changes the first particle into the new one with physics caluclated
    # before being sent to this function
    p1.v = v
    p1.type = flag  # sets new particle type

    # moves the other particles extremely far away with no velocity and with
    # type none, so it doesn't affect the rest of the particles
    # I did this to prevent having to remove the particles from
    # the lists.
    p2.type = None
    p2.v[0] = 0
    p2.v[1] = 0
    p2.x[0] = (np.random.random() + 500)*1000
    p2.x[0] = (np.random.random() + 500)*1000


def collisionDetection(pX, p, reaction):
    '''
    Detects all of the particles that are close enough to be considered
    colliding
    '''
    # makes an array with all of the distances between every combination
    # of particles.  Essentially finds all distances without using a loop
    dist = d.cdist(pX, pX)

    # defines index variables
    i = 0
    j = 0

    for i in range(len(p)):
        for j in range(i):
            # tests if each combination of the particles is close enough
            # to be colliding.  Uses particle diameters to determine
            # how close they can get (collideDist)
            collideDist = p[i].diameter/2 + p[j].diameter/2
            if dist[i, j] < collideDist:
                # if there is a collision, and no reaction is being
                # simulated, than a simple collision detection is done
                if not reaction:
                    collisionCalculaiton(p[i], p[j], dist[i, j], collideDist)
                # if there is a reaction between a red and blue particle,
                # it is carried out here
                elif (p[i].type == 'blue' and p[j].type == 'red') or        \
                     (p[i].type == 'red' and p[j].type == 'blue'):
                    flag = 'green'  # sets type of product of reaction
                    collisionReaction(p[i], p[j], i, j, p, pX, flag)
                # same for a reaction between blue and green particles
                elif (p[i].type == 'blue' and p[j].type == 'green') or      \
                     (p[i].type == 'green' and p[j].type == 'blue'):
                    flag = 'black'  # sets type of product of reaction
                    collisionReaction(p[i], p[j], i, j, p, pX, flag)
                # for all other combinations of particles, even during
                # a reaction, a normal collision is calculated
                else:
                    collisionCalculaiton(p[i], p[j], dist[i, j], collideDist)


def collisionReaction(p1, p2, i, j, p, pX, flag):
    '''
    This function caluclates the velocity of the product of two partcles
    chemically reacting with each other.  The flag variable is the end result
    of the reaction
    '''
    # particles masses
    m1 = p1.m
    m2 = p2.m

    # particle velocities
    v1 = p1.v
    v2 = p2.v

    # calculated velocity (assuming fully inelastic, this is the velocity at)
    # the center of mass
    vCMx = (m1*v1[0]+m2*v2[0])/(m1+m2)
    vCMy = (m1*v1[1]+m2*v2[1])/(m1+m2)

    # combines the two velocities into one array
    vCM = np.array((vCMx, vCMy))

    # sends data to function which carries out chemical reaction
    newParticle(p1, p2, vCM, flag)


def collisionCalculaiton(p1, p2, dist, collideDist):
    '''
    This does all of the calculations needed for an elastic collision between
    two non-reacting particles.
    '''

    # checks to see if the particles are close enough to collide
    delx = p1.x - p2.x
    distance = np.sqrt(np.sum(delx**2))
    # if there are in the exact same locatin, one is moved just slightly out
    # of the way
    if distance < p1.diameter/2 + p2.diameter/2:
        if distance == 0:
            p1.x -= 0.00001

        # orignally, when using the wikipedia equaiton, I missread the equaiton
        # and made the v-prime the new v, when it should have been deltaVself.
        # one of the sources showed me this mistake

        # this offset is added if the two particles are actually clipping into
        # each other
        offset = dist - (p1.diameter/2 + p2.diameter/2)
        p1.positionIncrease((-delx/distance)*offset/2)
        p2.positionIncrease((delx/distance)*offset/2)

        # elastic collision calculations
        total_m = p1.m + p2.m
        # changes in velocity
        delv1 = -2*p2.m/total_m*np.inner(p1.v-p2.v, p1.x-p2.x) /   \
            np.sum((p1.x-p2.x)**2)*(p1.x-p2.x)
        delv2 = -2*p1.m/total_m*np.inner(p2.v-p1.v, p2.x-p1.x) /    \
            np.sum((p2.x-p1.x)**2)*(p2.x-p1.x)
        # changes in velocity are added to each particle
        p1.velocityIncrease(delv1)
        p2.velocityIncrease(delv2)


class rectangleObject:
    '''
    This class simply keeps the data for the rectange that represents the
    container the particles are in.  Allows me to simulate the moving
    piston
    '''
    def __init__(self, dimentions):
        xDim, yDim = dimentions
        self.xDim = xDim
        self.yDim = yDim
        self.yDimLower = yDim


class particle:
    '''
    This is the main class for every particle created.
    '''
    def __init__(self, v, x, step, type):
        self.x = np.array(x)  # stores particle's current position
        self.v = np.array(v)  # current velocity
        self.diameter = 1     # diameter
        self.m = 1            # mass

        # the particle's type determines what it reacts with and what
        # color it shows up as.  Is a string of the color name ('red')
        self.type = type

    def move(self, i, p, pX, delt):
        '''
        moves the particle to the next location based on its velocity and
        the ammount of time each step reprents (delt)
        '''
        self.x += self.v*delt

    def wallCollision(self):
        '''
        calculates if the particle is hitting a wall and bounces it off of
        it if it is
        '''
        # gets dimentions of the tank it is located in
        xDim = rectOb.xDim
        yDim = rectOb.yDim  #changes when piston moves
        yDimLower = rectOb.yDimLower

        # renames particle location
        xpos = self.x[0]
        ypos = self.x[1]

        r = self.diameter / 2  # radius

        # checks if the particle is outside of each wall.  If it is, it
        # switches its velocity so it reflects off the wall.
        if xpos + r >= xDim:
            # this check is done for each one to see if the particle is
            # already moving away from the wall.  Without this check, the
            # partcle's velocity can keep swtiching forever, making the
            # stuck in the wall
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
        '''
        sets particle velocity to given one
        '''
        self.v = v

    def positionIncrease(self, dx):
        '''
        increse position with given increase
        '''
        self.x = self.x + dx

    def velocityIncrease(self, dv):
        '''
        increases velocity given velocity increase
        '''
        self.v = self.v + dv


class trackedParticle(particle):
    '''
    subclass of particle, used to add a larger tracer particle that
    leaves behind a trail, showing its movement throughout the simulaiton
    '''

    def __init__(self, x, v, step, type, steps):
        particle.__init__(self, v, x, step, type)  # normal particle physics
        # records all of its positions, allowing a line of its path to be
        # printed
        self.xRecord = np.empty((steps, 2), dtype=float)
        self.index = 0
        # has a larger diameter so it is easier to see
        self.diameter = 1.6

    def move(self, i, p, pX, delt):
        '''
        in additon to moving normally, it also records all of its past
        locations, allowing its path to be printed on the graph
        '''
        self.x += self.v*delt
        self.xRecord[self.index] = self.x  # record of all positions
        self.index += 1  # increases index by one, so it rights new positions


def runsim(mode='double', reaction=True, seed=0, nParticles=20, dimentions=(10, 10), piston=True, steps=500, track=True, baseline=None, displacement = 0.01):
    '''
    This program runs the full simulation, giving a graphical output
    '''
    # makes the rectangle object global, so the slider bar can adjust it
    global rectOb
    rectOb = rectangleObject(dimentions)

    xDim, yDim = dimentions

    # sets random seed for consistant results
    np.random.seed(seed)

    # if two particles are being used, this breaks the nParticles up so that
    # the amount of each type is set
    if mode == 'double':
        nParticlesBlue = int(nParticles/2)
        nParticlesRed = int(nParticles-nParticlesBlue)
        nParticles = nParticlesBlue + nParticlesRed
    else:
        nParticlesBlue = nParticles-1

    step = 40

    # sets up the three main data carriers (all empty for now)
    # p is a list of all of the particle objects, so the are easily indexed
    p = [0]*nParticles
    # pX is all of the postions in an array, for easy calculations
    pX = np.zeros((nParticles, 2))
    # pT is all of the types in an array, for easy printing (used in graphing)
    pT = np.empty((nParticles), dtype=object)

    # sets up all of the blue particles
    for j in range(nParticlesBlue+1):
        # randomly decides each of their x and v values
        x = np.array((float(np.random.rand(1)),
                      float(np.random.rand(1))))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))
        # creates an isntance of the particle class
        p[j] = particle(v, x, step, 'blue')

        # records all of their positions
        pX[j, 0] = p[j].x[0]
        pX[j, 1] = p[j].x[1]

    # if two particles are being used, this adds the red ones
    if mode == 'double':
        for k in range(nParticlesRed):
            k += j
            x = np.array((float(np.random.rand(1))*8,
                          float(np.random.rand(1))*8))
            v = np.array((float(np.random.rand(1))*50,
                          float(np.random.rand(1))*50))
            p[k] = particle(v, x, step, 'red')
            pX[k, 0] = p[k].x[0]
            pX[k, 1] = p[k].x[1]

    # if a tracking particle is being used, it is added here
    if track:
        x = np.array((float(np.random.rand(1)+displacement),
                      float(np.random.rand(1)+displacement) ))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))

        # creates instance of trackedParticle
        tracked = trackedParticle(x, v, step, 'tracked', steps)

        # adds the instance, its positon, and its type to the correct locations
        p = np.append(p, tracked)
        pX = np.vstack([pX, tracked.x])
        pT = np.append(pT, tracked.type)
        # increase the total nParticles by one
        nParticles += 1

    # sets up the plots used for animation
    fig = plt.figure()
    ax = plt.axes()

    # defines all the points and lines being animated
    point1, = ax.plot([], [], 'bo', markersize=8)  # blue particles
    point2, = ax.plot([], [], 'ro', markersize=8)  # red
    point3, = ax.plot([], [], 'go', markersize=8)  # green
    point4, = ax.plot([], [], 'ko', markersize=8)  # black
    point5, = ax.plot([], [], 'yo', markersize=14)  # tracked particle
    linea,   = ax.plot([], [], 'orange')  # current tracked particle path
    lineb,  = ax.plot([], [], 'c:')  # previous tracked particle

    # sets plot dimentions
    ax.set_xlim(-xDim-1, xDim+1)
    ax.set_ylim(-yDim-1, yDim+1)
    plt.title('Simulating particles moving in a box')
    plt.legend((linea, lineb), ('Current Path (with added \
starting displacement)', 'Baseline Path'))
    # defines the rectangle representing the box they are in
    rectCentX, rectCentY = -(xDim), -(yDim)
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((rectCentX, rectCentY), 2*xDim, 2*yDim,
                          fill=False))

    # if the piston is desired, a new line, representing it is printed
    if piston:
        x = (-xDim, xDim)
        y = (yDim, yDim)
        [line] = ax.plot(x, y)  # draws the line
        # Define an axes area and draw a slider in it
        # This slider allows the height of the box to be adjusted
        # during the simulation
        amp_slider_ax = fig.add_axes([0.25, 0.02, 0.55, 0.03])
        amp_slider = Slider(amp_slider_ax, 'Height', -yDim, yDim,
                            valinit=rectOb.yDim)

        # Define an action for modifying the line when any slider value changes
        def sliders_on_changed(val):
            line.set_ydata(amp_slider.val)
            fig.canvas.draw_idle()
            rectOb.yDim = val  # moves the top of the box to slider position
        amp_slider.on_changed(sliders_on_changed)

    # change in time per update
    dt = 0.01
    abcde = 0  # counter used to limit the number of updates
    while abcde < steps:  # limits updates to givin number of steps
        for types in range(len(p)):
            # makes array of each particle's type for indexing
            pT[types] = p[types].type

        # find each particle thats blue
        blueIndex = pT == 'blue'

        # graphs each blue particle, using index calculated above
        xvalBlue = pX[blueIndex, 0]
        yvalBlue = pX[blueIndex, 1]
        point1.set_data(xvalBlue, yvalBlue)

        # if multiple types of particles, they are added as well
        if mode == 'double':

            # adds red particles
            redIndex = pT == 'red'
            xvalRed = pX[redIndex, 0]
            yvalRed = pX[redIndex, 1]
            point2.set_data(xvalRed, yvalRed)

            if reaction:
                # graphs green particles if a reaction is taking place
                greenIndex = pT == 'green'
                xvalGreen = pX[greenIndex, 0]
                yvalGreen = pX[greenIndex, 1]
                point3.set_data(xvalGreen, yvalGreen)

                # black particles
                blackIndex = pT == 'black'
                xvalBlack = pX[blackIndex, 0]
                yvalBlack = pX[blackIndex, 1]
                point4.set_data(xvalBlack, yvalBlack)

        # if there is a tracking particle, it is graphed, in addition to
        # its path
        if track:
            xvalTrack = pX[-1, 0]
            yvalTrack = pX[-1, 1]
            point5.set_data(xvalTrack, yvalTrack)  # plots tracker

            recIndex = p[-1].index
            xRec = p[-1].xRecord[:recIndex, 0]
            yRec = p[-1].xRecord[:recIndex, 1]
            linea.set_data(xRec, yRec)  # plots path

            # if baseline data is given (data from a previos tracker particle)
            # it's path is tracked as well.  This is used for comparison
            # purposes
            if baseline != None:
                xbase, ybase = baseline
                lineb.set_data(xbase[:recIndex], ybase[:recIndex])

        # updates every particle
        for i in range(len(p)):
            p[i].move(i, p, pX, dt)  # moves each to new location
            p[i].wallCollision()  # check and calculate wall collision
            collisionDetection(pX, p, reaction)  # particle-particle collision
        for j in range(nParticles):
            # updates position information
            pX[j] = p[j].x
        plt.pause(.0001)  # waits a very short bit
        abcde += 1  # counts this step

def runsimCompare(mode='double', reaction=True, seed=0, nParticles=20, dimentions=(10, 10), piston=True, steps=500, track=True):
    '''
    does same simulation as runsim(), but does not give any graphical output
    Instead, it returns the path of the tracked particle to be graphed on
    by the runsim() funciton.  See runsim() for documention on how the
    funciton works
    '''
    global rectOb
    rectOb = rectangleObject(dimentions)

    xDim, yDim = dimentions
    np.random.seed(0)
    if mode == 'double':
        nParticlesBlue = int(nParticles/2)
        nParticlesRed = int(nParticles-nParticlesBlue)
        nParticles = nParticlesBlue + nParticlesRed
    else:
        nParticlesBlue = nParticles-1

    step = 40
    # v = (5, 5)
    p = [0]*nParticles
    pX = np.zeros((nParticles, 2))
    pT = np.empty((nParticles), dtype=object)

    for j in range(nParticlesBlue+1):
        x = np.array((float(np.random.rand(1)),
                      float(np.random.rand(1))))
        v = np.array((float(np.random.rand(1))*50,
                      float(np.random.rand(1))*50))
        # create an instance of the drunkard class
        p[j] = particle(v, x, step, 'blue')
        pX[j, 0] = p[j].x[0]
        pX[j, 1] = p[j].x[1]

    if mode == 'double':
        for k in range(nParticlesRed):
            k += j
            x = np.array((float(np.random.rand(1))*8,
                          float(np.random.rand(1))*8))
            v = np.array((float(np.random.rand(1))*50,
                          float(np.random.rand(1))*50))
            p[k] = particle(v, x, step, 'red')
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
        nParticles += 1


    dt = 0.01
    abcde = 0
    while abcde < steps:
        for types in range(len(p)):
            pT[types] = p[types].type

        for i in range(len(p)):
            p[i].move(i, p, pX, dt)
            p[i].wallCollision()
            collisionDetection(pX, p, reaction)
        for j in range(nParticles):
            pX[j] = p[j].x

        abcde += 1


    return p[-1]

def presentation(mode='double', reaction=True, seed=0, nParticles=20, dimentions=(10, 10), piston=True, steps=500, track=True, displacement=0.0000001, compare=True):
    '''
    This function simulates a group of particles in a box, including collision.
    It is intended to be both a visual representation of particle movement
    as well as showing how chaotic this type of system can be.

    NOTE: Running the simulation with too many particles makes it very slow
    NOTE: Decreasing the box height too quickly can cause the particles to be
    temporarily stuck outside the box.
    NOTE: Can take some time for the animation to show if up comparing, since
    first simulation must be completed

    PARAMETERS:
        mode: either 'singe' or 'double'.  Single adds only one type of
            particle while doubles adds two types
        reaction: either true or false, tell the program whether a reaction
            between the particles takes place.  The reaction is:

                red + blue ------> green
                green + blue ----> black

        seed: sets the seed for the random particle position and velocity
            starting values
        nParticles: Integer setting the total number of particle
        dimentions: Tuple containing the hight and width of the containing box
        piston: True or False, determines of the piston is added.  Adding the
            piston adds both a line representing the new box height, and
            a slider.  One can click on the slider to move the piston up or
            down, allowing the height of the box to be adjustable during
            the simulation
        steps: An integer setting the number of update interations carried out
        track: True or False, tells the program whether or not to include the
            tracking particle.  The tracking particle is of a different color
            and leaves its path behind, allowing one to better visualize the
            particle movement
        displacement: If compare is True, sets the amount of distance the
            tracking particle starting position differs between the first
            simulation and the second.  Is a float
        compare: True or False.  If true, the simulation is run twice. The
            first time, the simulation gives no graph. Instead, the path of
            the tracked particle from the first simulation is printed onto the
            graph of the second.  This way, the path of the tracker particle
            in the second simulation can becompared to that of the first.

            The first acts as a baseline.  When the second is run, the
            staring location of the tracker particle in the second simulation
            is moved by the amount specified as the displacement.

    OUTPUT:
        Program gives a simulation of the particles given the environment
        specified above.
    '''
    if mode != 'single' and mode != 'double':
        raise Exception('mode must be single or double')
    if reaction != True and reaction != False:
        raise Exception('reaction can only be true or false')
    if type(seed) != int and type(seed) != float:
        raise Exception('seed must be a number')
    if type(nParticles) != int:
        raise Exception('nParticles must be an integer')
    if type(dimentions) != tuple:
        raise Exception('dimentions must be in a tuple')
    if piston != True and piston != False:
        raise Exception('piston is either True or False')
    if type(steps) != int:
        raise Exception('steps must be an integer')
    if track != True and track != False:
        raise Exception('track must be True or False')
    if type(displacement) != float and type(displacement) != int:
        raise Exception('displacement must be a number')
    if compare != True and compare != False:
        raise Exception('compare is either True or False')

    if compare:
        # if there is a comparison to be made, the path of the tracking
        # particle from the first simulation is found here.
        pX = runsimCompare(mode, reaction, seed, nParticles, dimentions, piston, steps, track)
        xplot = pX.xRecord[:, 0]
        yplot = pX.xRecord[:, 1]

        # the results from the first simulation is passed to the second one
        # so that it can be printed onto the animation
        baseline = (xplot, yplot)
    else:
        baseline = None  # if no comparison, baseline set to none to be ignored
        displacement = 0 # keeps tracker at inital start position

    # runs graphical simulation
    runsim(mode, reaction, seed, nParticles, dimentions, piston, steps, track, baseline, displacement)
