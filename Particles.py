import numpy as np

damp = 0.01
# define particle class
class Particle():
    def __init__(self, env, X, V, A, IntForce, ExtForce, radius, H, mass, density, volume, pressure, color, RestPosition, neighbors, Plastic_epsilon):
        self.env = env
        self.X = X
        self.V = V
        self.A = A
        self.IntForce = IntForce
        self.ExtForce = ExtForce
        self.r = radius
        self.H = H
        self.m = mass
        self.rho = density
        self.vol = volume
        self.p = pressure
        self.color = color
        self.ResPos = RestPosition
        self.neighbors = neighbors
        self.closestParticle = None
        self.PrevPos = None
        self.NumOfNeigh = None
        self.MomentMatrix = None
        self.Fixed = False
        self.u = None
        self.B = None
        self.Du = None
        self.Jacobian = None
        self.Plastic_epsilon = Plastic_epsilon
        self.Strain_epsilon = None
        self.Stress_Sigma = None
        self.Strain_epsilonP = None
        self.Fe = None
        self.Fv = None
        self.F_V_E = None
        self.F_V_E_invA = None
        self.F_press = None
        self.F_visc = None


    def IsFixed(self):
        self.Fixed = True

    def UpdateAcceleration(self):
        if self.Fixed != True:
            self.A = (self.IntForce + self.ExtForce) / self.rho  # denominator was mass instead of density
        else:
            self.A = np.array([0, 0, 0], dtype=float)

    def addAcceleration(self, acc):
        self.A += acc

    def addVelocity(self, vel):
        self.V += vel

    def addPosition(self, pos):
        self.X += pos

    def attract(self, particle):
        if particle not in self.neighbors:
            r = self.X - particle.X
            # self.A += 6.67408e-11 * particle.m / r ** 2
            particle.ExtForce += particle.m * (0.1 * particle.m / np.linalg.norm(r)) * r

    def UpdateVelocity(self):
        self.V += self.A * self.env.dt
        self.V *= (1 - damp)

    def UpdatePosition(self):
        self.X += self.V * self.env.dt #- 0.5 * self.A * self.env.dt ** 2