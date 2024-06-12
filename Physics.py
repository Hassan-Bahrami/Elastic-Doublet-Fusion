import Particles as pf
import matplotlib.cm as cm
import numpy as np
from scipy.spatial import KDTree



import cProfile, pstats, io


def profile(fnc):
    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner

np.set_printoptions(suppress=True, precision=7)


# Solver parameters
damp = 0.01

# Smoothing kernels defined in Müller and their gradients
# H = 50
# Poly6 = 315 / (64 * np.pi * (H ** 9))
# SpikyGrad = -45 / (np.pi * (H ** 6))
# ViscLapla = 45 / (np.pi * (H ** 6))
# k_v = 1

# Model parameters
e = 0.5
nu = 0.49
kv = 1
s = 49.745623
ro_uni = 10
Gamma1 = 0.0001
Gamma2 = 1.5 # changed from 0.1 to 1.5
G = 3

class Environment():
    def __init__(self, GRAVITY, dt):
        self.GRAVITY = GRAVITY
        self.dt = dt
        self.Sphere1 = []
        self.Sphere2 = []
        self.particles = []
        self.mesh = None
        self.E = None
        self.Nu = None
    # @profile
    def update_physics(self):
        self.E = e
        self.Nu = nu
        for p1 in self.particles:
            self.Reset(p1)

        # tree = KDTree(self.update_PointCloud())
        # for p1 in self.particles:
        #     self.update_neighbor(p1, tree)

        for p1 in self.particles:
            self.update_density(p1)
            self.update_volume(p1)

        for p1 in self.particles:
            self.Inv_MomentMatrix(p1)
        center_of_mass1 = self.center_of_mass(self.Sphere1)
        center_of_mass2 = self.center_of_mass(self.Sphere2)
        center_difference = np.linalg.norm(center_of_mass2 - center_of_mass1)
        if center_difference > 1.5:
            for p1 in self.Sphere1:
                self.compute_external_Force(p1, center_of_mass2, G=3, avg_wght=1.1, mu=0.1, A=5)
            for p2 in self.Sphere2:
                self.compute_external_Force(p2, center_of_mass1, G=3, avg_wght=1.1, mu=0.1, A=5)

        for p1 in self.Sphere1:
            for p2 in self.Sphere2:
                self.elasticCollision(p1, p2)



        for i, p1 in enumerate(self.particles):
            self.Compute_internalForce(p1)

        for i, p1 in enumerate(self.particles):
            p1.UpdateAcceleration()

        for i, p1 in enumerate(self.particles):
            p1.UpdateVelocity()

        for i, p1 in enumerate(self.particles):
            p1.UpdatePosition()


    def addParticle2Sphere1(self, p):
        self.Sphere1.append(p)
        p.addAcceleration(self.GRAVITY)

    def addParticle2Sphere2(self, p):
        self.Sphere2.append(p)
        p.addAcceleration(self.GRAVITY)

    def addParticle(self, p):
        self.particles.append(p)
        p.addAcceleration(self.GRAVITY)

    # def bounce(self, p):
    #     i = 0
    #     for x in p.X[0]:
    #         if x > self.DIM[i] - p.r:
    #             dist = p.r - (self.DIM[i] - x)
    #             p.addPosition(-dist)
    #             tmp = np.zeros(np.size(p.V))
    #             tmp[i] = -2 * p.V[0][i]
    #             p.addVelocity(tmp * 0.8)
    #         elif x < p.r:
    #             dist = p.r - x
    #             p.addPosition(dist)
    #             tmp = np.zeros(np.size(p.X))
    #             tmp[i] = -2 * p.V[0][i]
    #             p.addVelocity(tmp * 0.8)
    #         i += 1

    def elasticCollision(self, p1, p2):
        dX = p1.X-p2.X
        dist = np.sqrt(np.sum(dX**2))
        if dist < p1.r + p2.r:
            offset = dist-(p1.r+p2.r)
            p1.addPosition((-dX/dist)*offset/2)
            p2.addPosition((dX/dist)*offset/2)
            total_mass = p1.m+p2.m
            dv1 = -2*p2.m/total_mass*np.inner(p1.V-p2.V,p1.X-p2.X)/np.sum((p1.X-p2.X)**2)*(p1.X-p2.X)
            dv2 = -2*p1.m/total_mass*np.inner(p2.V-p1.V,p2.X-p1.X)/np.sum((p2.X-p1.X)**2)*(p2.X-p1.X)
            p1.addVelocity(dv1)
            p2.addVelocity(dv2)

    def plasticCollision(self):
        pass

    # computing kernel radius of each particle
    def finding_Kernel_radius_Sphere1(self, sphere, k):
        # sphere1 = point cloud
        # k = number of nearest neighbours
        tree = KDTree(sphere)
        for i, p1 in enumerate(self.Sphere1):
            dis, idx = tree.query(p1.X, k)
            p1.H = 1.1 * np.max(dis)

    def finding_Kernel_radius_Sphere2(self, sphere, k):
        # sphere1 = point cloud
        # k = number of nearest neighbours
        tree = KDTree(sphere)
        for i, p1 in enumerate(self.Sphere2):
            dis, idx = tree.query(p1.X, k)
            p1.H = 1.1 * np.max(dis)

    # find the nearest neighbor of each particle
    def NeighborsOfSphere1(self, sphere):
        tree = KDTree(sphere)
        for i, p1 in enumerate(self.Sphere1):
            p1.neighbors = []
            idx = tree.query_ball_point(p1.X, p1.H)
            for j in idx:
                p2 = self.Sphere1[j]
                if p1 != p2:
                    p1.neighbors.append(p2)

    def NeighborsOfSphere2(self, sphere):
        tree = KDTree(sphere)
        for i, p1 in enumerate(self.Sphere2):
            p1.neighbors = []
            idx = tree.query_ball_point(p1.X, p1.H)
            for j in idx:
                p2 = self.Sphere2[j]
                if p1 != p2:
                    p1.neighbors.append(p2)

    def find_neighbors(self, sphere):
        positions = [particle.X for particle in self.Sphere1]
        tree = KDTree(positions)

        for i, particle in enumerate(self.Sphere1):
            particle.neighbors = []
            indices = tree.query_ball_point(particle.X, particle.H)
            particle.neighbors = [j for j in indices if j != i]

    def find_neighbors_for_spheres(self):
        self.find_neighbors(self.Sphere1)
        self.find_neighbors(self.Sphere2)

    def colosest_particle_in_spehriod2(self, sphere):
        tree = KDTree(sphere)
        for p1 in self.Sphere1:
            dist, ind = tree.query(p1.X, k=1)
            cord1 = tree.data[ind]
            for p2 in self.Sphere2:
                if np.isclose(p2.X, cord1, atol=1e-3).all():
                    p1.closestParticle = p2
                    break


    def update_neighbor(self, p1, tree):
        p1.neighbors = []
        idx = tree.query_ball_point(p1.X, p1.H)
        for j in idx:
            p2 = self.particles[j]
            if p1 != p2:
                p1.neighbors.append(p2)

    def update_PointCloud(self):
        points = np.zeros((len(self.particles), 3))
        for i, p1 in enumerate(self.particles):
            points[i] = p1.X
        return points

    # compute mass of each particle
    def Compute_ParticleMassOfSphere1(self, s, rho_uni):
        # s = scaling factor (same for all particles)
        # rho = material density
        for i, p1 in enumerate(self.Sphere1):
            p1.m = s * ((p1.H/3) ** 3) * rho_uni

    def Compute_ParticleMassOfSphere2(self, s, rho_uni):
        # s = scaling factor (same for all particles)
        # rho = material density
        for i, p1 in enumerate(self.Sphere2):
            p1.m = s * ((p1.H/3) ** 3) * rho_uni


    # apply color map
    def ColorMapOfSphere1(self):
        viridis = cm.get_cmap('plasma', 20)
        for i, p1 in enumerate(self.Sphere1):
            p1.color = np.asarray(viridis(p1.m))

    # apply color map
    def ColorMapOfSphere2(self):
        viridis = cm.get_cmap('plasma', 20)
        for i, p1 in enumerate(self.Sphere2):
            p1.color = np.asarray(viridis(p1.m))

    # weighting function
    def Wij(self, r, h):
        norm_term = 4 * np.pi * h*h/9
        if r >= h:
            return 0
        return (315 / (64 * np.pi * (h ** 9))) * ((h*h - r*r) ** 3) * norm_term

    # compute density of each particle
    def Compute_densityOfSphere1(self):
        for i, p1 in enumerate(self.Sphere1):
            p1.rho = 0
            for p2 in p1.neighbors:
                if p1 == p2:
                    continue
                dX = p2.X - p1.X
                #r = np.sqrt(np.sum(dX ** 2))
                r = np.linalg.norm(dX)
                p1.rho += p2.m * self.Wij(r, p1.H)

    # compute density of each particle
    def Compute_densityOfSphere2(self):
        for i, p1 in enumerate(self.Sphere2):
            p1.rho = 0
            for p2 in p1.neighbors:
                if p1 == p2:
                    continue
                dX = p2.X - p1.X
                #r = np.sqrt(np.sum(dX ** 2))
                r = np.linalg.norm(dX)
                p1.rho += p2.m * self.Wij(r, p1.H)

    def update_density(self, p1):
        p1.rho = 0
        for p2 in p1.neighbors:
            if p1 == p2:
                continue
            dX = p2.X - p1.X
            # r = np.sqrt(np.sum(dX ** 2))
            r = np.linalg.norm(dX)
            p1.rho += p2.m * self.Wij(r, p1.H)

    # Compute particle volumes
    def Compute_volumeOfSphere1(self):
        for i, p1 in enumerate(self.Sphere1):
            p1.vol = p1.m/p1.rho

    # Compute particle volumes
    def Compute_volumeOfSphere2(self):
        for i, p1 in enumerate(self.Sphere2):
            p1.vol = p1.m/p1.rho

    def update_volume(self, p1):
        p1.vol = p1.m / p1.rho

    def Inv_MomentMatrix(self, p1):
        # for i, p1 in enumerate(self.particles):
        A = np.zeros((3, 3))
        for j, p2 in enumerate(p1.neighbors):
            x_ij = p2.ResPos - p1.ResPos
            #r = np.sqrt(np.sum(x_ij ** 2))
            r = np.linalg.norm(x_ij)
            w_ij = self.Wij(r, p1.H)
            A_sum = np.outer(x_ij, x_ij)
            A_sum *= w_ij
            A += A_sum
            p1.B = A
        if np.linalg.det(p1.B) != 0:
            inv_A = np.linalg.inv(p1.B)
        else:
            inv_A = np.identity(3)
        p1.MomentMatrix = inv_A

    def update_Inv_MomentMatrix(self, p1, p2):
        A = np.zeros((3, 3))
        A_sum = np.zeros((3, 3))
        x_ij = p2.X - p1.X
        # r = np.sqrt(np.sum(x_ij ** 2))
        r = np.linalg.norm(x_ij)
        w_ij = self.Wij(r, p1.H)
        A_sum = np.outer(x_ij, x_ij)
        A_sum *= w_ij
        p1.B += A_sum
        p1.MomentMatrix = np.linalg.inv(p1.B)

    def Reset(self, p1):
        p1.IntForce = np.zeros((3,))
        p1.ExtForce = np.zeros((3,))
        p1.A = np.zeros((3,))
        p1.u = p1.X - p1.ResPos
        p1.PrevPos = p1.X
        # p1.Plastic_epsilon = np.zeros((3, 3))

    def Compute_PlasticStrain(self, epsilon, epsilonP, Gamma1, Gamma2):
        epsilon_e = epsilon - epsilonP
        epsilon_prime = epsilon_e - ((np.trace(epsilon_e)/3) * np.identity(3))
        N1 = np.linalg.norm(epsilon_prime)
        if N1 > Gamma1:
            deltaEpsilonP = ((N1 - Gamma1)/N1) * epsilon_prime
            epsilonP1 = epsilonP + deltaEpsilonP
            N2 = np.linalg.norm(epsilonP1)
            epsilonP = epsilonP1 * np.minimum(1, Gamma2/N2)
        return epsilonP

    def Compute_Stress(self, epsilon, epsilonP, Mu, Lambda):
        sigma = np.zeros((3, 3))
        epsilon1 = epsilon - epsilonP
        Tr = np.trace(epsilon1)
        for i in range(3):
            for j in range(i+1):
                sigma[i][j] = epsilon1[i][j] * 2.0 * Mu
            sigma[i][i] += Lambda * Tr
        # Symmetrize sigma
        sigma[0][1] = sigma[1][0]
        sigma[0][2] = sigma[2][0]
        sigma[1][2] = sigma[2][1]
        epsilon = epsilon1
        return epsilon, sigma

    def Compute_Fe(self, J, sigma, v):
        Fe = (J @ sigma) * 2 * (-v)
        return Fe

    def Compute_Fv(self, J, Lambda, kv, v):
        det_J = np.linalg.det(J)
        det_J_derivative = np.zeros((3, 3))
        det_J_derivative[0][:] = np.cross(J[1, :], J[2, :])
        det_J_derivative[1][:] = np.cross(J[2, :], J[0, :])
        det_J_derivative[2][:] = np.cross(J[0, :], J[1, :])
        Fv = det_J_derivative * (det_J - 1) * Lambda * kv * (-v)
        return Fv

    def Set_LambdaMu(self, E, Nu):
        Lambda = Nu * E / ((1 + Nu) * (1 - 2 * Nu))
        Mu = E / 2 * (1 + Nu)
        return Mu, Lambda

    def Compute_internalForce(self, p1):
        B = np.zeros((3, 3))
        Du = np.zeros((3, 3))
        J = np.zeros((3, 3))
        epsilon = np.zeros((3, 3))
        sigma = np.zeros((3, 3))
        B_sum = np.zeros((3, 3))
        J_tr = np.zeros((3, 3))
        # for i, p1 in enumerate(self.particles):
        for j, p2 in enumerate(p1.neighbors):
            x_ij = p2.ResPos - p1.ResPos
            r = np.sqrt(np.sum(x_ij ** 2))
            u_ij = p2.u - p1.u #u_j - u_i
            B_sum = np.outer(u_ij, x_ij)
            B_sum *= self.Wij(r, p1.H)
            B += B_sum
        B = np.transpose(B)
        p1.B = B
        Du = p1.MomentMatrix @ B
        Du = np.transpose(Du)
        p1.Du = Du
        I = np.identity(3)
        J = Du + I
        p1.Jacobian = J
        J_tr = np.transpose(J)
        epsilon = (J @ J_tr) - I
        p1.Strain_epsilon = epsilon
        # computePlasticStrain
        if Gamma1 > 0:
            epsilonP = self.Compute_PlasticStrain(epsilon=p1.Strain_epsilon, epsilonP=p1.Plastic_epsilon, Gamma1=Gamma1, Gamma2=Gamma2)
            p1.Plastic_epsilon = epsilonP
        # Compute Mu and Lambda
        Mu, Lambda = self.Set_LambdaMu(E=self.E , Nu=self.Nu)
        # computeStress
        epsilon, sigma = self.Compute_Stress(epsilon, epsilonP=p1.Plastic_epsilon, Mu=Mu, Lambda=Lambda)
        p1.Stress_Sigma = sigma
        p1.Strain_epsilon = epsilon
        # Compute Fe
        Fe = self.Compute_Fe(J, p1.Stress_Sigma, v=p1.vol)
        p1.Fe = Fe
        # Compute Fv
        Fv = self.Compute_Fv(J, Lambda=Lambda, kv=kv, v=p1.vol)
        p1.Fv = Fv
        F_V_E = Fv + Fe
        p1.F_V_E = F_V_E
        F_V_E_invA = F_V_E @ p1.MomentMatrix
        p1.F_V_E_invA = F_V_E_invA
        for j, p2 in enumerate(p1.neighbors):
            force_to_j = np.zeros((3,3))
            x_ij = p2.ResPos - p1.ResPos
            r = np.sqrt(np.sum(x_ij ** 2))
            force_to_j = F_V_E_invA @ (x_ij * self.Wij(r, p1.H))
            if p2.Fixed != True:
                p2.IntForce += force_to_j
            p1.IntForce -= force_to_j

    def extract_vertsOfSphere1(self):
        vertsOfSphere1 =[]
        for p in self.particles[0:int((len(self.particles)/2))]:
            vertsOfSphere1.append(p.X)
        return vertsOfSphere1

    def extract_vertsOfSphere2(self):
        vertsOfSphere2 =[]
        for p in self.particles[int((len(self.particles)/2)):len(self.particles)+1]:
            vertsOfSphere2.append(p.X)
        return vertsOfSphere2

    def extract_colorsOfSphere1(self):
        colorsOfSphere1 =[]
        for p in self.Sphere1:
            colorsOfSphere1.append(p.color)
        return colorsOfSphere1

    def extract_colorsOfSphere2(self):
        colorsOfSphere2 =[]
        for p in self.Sphere2:
            colorsOfSphere2.append(p.color)
        return colorsOfSphere2

    def attract(self, p1, p2):
        G = 5
        mu = 1.0016
        A = 25
        if p1 not in p2.neighbors:
            r = p2.X - p1.X
            p1.ExtForce = np.array([0, 0, 0], dtype=float)
            # self.A += 6.67408e-11 * particle.m / r ** 2
            # p1.ExtForce += p1.m * (0.1 * p2.m / np.linalg.norm(r)) * r
            Gravity_force = (G * p1.m * p2.m / np.linalg.norm(r)**2) * r / np.linalg.norm(r)
            Viscosity_force = mu * p1.V * A * Gravity_force / np.linalg.norm(Gravity_force)
            p1.ExtForce += Gravity_force + Viscosity_force

    def compute_gravity_force(self, p1, G, avg_wght, r):
        gravity_force = G * p1.m * avg_wght / r**2
        return gravity_force

    def compute_external_Force(self, p1, center_of_mass, G, avg_wght, mu, A):
        # p2 = p1.closestParticle
        r = center_of_mass - p1.X
        gravity_force = G * p1.m * avg_wght
        viscosity_force = mu * p1.V * A
        total = gravity_force - viscosity_force
        p1.ExtForce += total * (r/np.linalg.norm(r))
        # p2.ExtForce += -total * (r/np.linalg.norm(r))

    def reset_RestPosition(self):
        for p in self.particles:
            p.ResPos = p.X

    def collision(self):
        # self.reset_RestPosition()
        tree = KDTree(self.update_PointCloud())
        for p1 in self.particles:
            self.update_neighbor(p1, tree)
        for p1 in self.particles:
            self.Inv_MomentMatrix(p1)

    def surface_extraction(self,c, rez):
        tree = KDTree(self.update_PointCloud())
        X, Y, Z = np.indices((20, 20, 20)) / rez
        surface = -c
        for x in X:
            for y in Y:
                for z in Z:
                    for p1 in self.particles:
                        dX = p1.X - [x, y, z]
                        # r = np.sqrt(np.sum(dX ** 2))
                        r = np.linalg.norm(dX)
                        if r <= p1.H:
                            surface[x, y, z] += p1.v * self.Wij(r, p1.H)

    def average_weight(self, sphere):
        mass_sum = 0
        counter = 0
        for p1 in sphere:
            mass_sum += p1.m
            counter += 1
        avg_wght = mass_sum / counter
        return avg_wght

    # finding the center mass of particle system
    def center_of_mass(self, sphere):
        Weighted_position = np.zeros((3,))
        total_mass = 0
        # counter = 0
        for p in sphere:
            Weighted_position += p.m * p.X
            total_mass += p.m
            # counter += 1
        center_of_mass = Weighted_position / total_mass
        # avg_mass = total_mass / counter
        return center_of_mass

    def compute_density_pressure(self, p1):
        p1.rho = 0

        for p2 in p1.neighbors:
            xij = p2.X - p1.X
            r = np.linalg.norm(xij)
            w_ij = self.Wij(r, p1.H)
            p1.rho += p2.m * w_ij

