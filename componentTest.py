import numpy as np

def distance(particle1,particle2):
    dist = np.sqrt( (particle1.xpos - particle2.xpos)**2 + (particle1.ypos - particle2.ypos))

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

    def collisionDetection(self, i,p):
        dist = np.zeros(len(p))
        for j in range(len(p)):
            if j != i:
                dist[j] = distance(p[i],p[j])
                if dist[j] < self.diameter:
                    r1 = np.array((self.xpos,self.ypos))
                    r2 = np.array((p[j].xpos,p[j].ypos))

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

        #Velocity in feet per second

        if self.xpos >= 10:
            if self.ypos >= 10:
                thetaRange = (np.pi, 3*np.pi/2)

            elif self.ypos <= -10:
                thetaRange = (np.pi/2 , np.pi)
            else:
                thetaRange = (np.pi/2 , 3*np.pi/2)
        elif self.ypos >= 10:
            thetaRange = (np.pi , 2*np.pi)

        elif self.xpos <= -10:
            if self.ypos >= 10:
                thetaRange = (3*np.pi/2, 2*np.pi)
            elif self.ypos <= -10:
                thetaRange = (0,np.pi/2)
            else:
                thetaRange = (3*np.pi/2 , 5*np.pi/2)
        elif self.ypos <= -10:
            thetaRange = (0 , np.pi)


        if abs(self.xpos) >= 10 or abs(self.ypos) >= 10:
            thetaHigh , thetaLow = thetaRange
            #the following comes from plugging in generic high and low values
            #into point-slop formula for a line
            self.theta = (thetaHigh - thetaLow) * np.random.random() + thetaLow
        self.xpos = self.v*np.cos(self.theta) + self.xpos
        self.ypos = self.v*np.sin(self.theta) + self.ypos

        self.batteryLife -= 1

        self.collisionDetection(i,p)
    def speed(self, v):
        self.v = v

particle1 = robovac(40,1)
particle2 = robovac(40,5)

p = [particle1, particle2]
dist = np.zeros(2)
dist[0] = distance(p[0],p[1])
print(dist)
