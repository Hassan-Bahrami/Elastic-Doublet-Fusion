import Physics as ph
import Particles as pt
import numpy as np
from spherical_sampling import fibonacci

# Environment parameter
GRAVITY = np.asarray([0, 0, 0])
dt = 0.01
env = ph.Environment(GRAVITY, dt)

s3 = np.zeros((1000, 3))
c = 0
s4 = []
for l in np.linspace(-1, 1, 5):
    for m in np.linspace(-1, 1, 5):
        for n in np.linspace(-1, 1, 5):
            s3[c] = np.array([l, m, n])
            if np.linalg.norm(s3[c]) < 1:
                s4.append(s3[c])
                c += 1
s4 = np.array(s4)
N = 30
samples = fibonacci(N, co_ords='cart')
merged_s4 = np.concatenate((s4, samples), axis=0)

s7 = np.zeros((len(merged_s4), 3))
s8 = np.zeros((len(merged_s4), 3))
for i in range(len(s7)):
    s7[i][0] = merged_s4[i][0]
    s7[i][1] = merged_s4[i][1] + 1.25
    s7[i][2] = merged_s4[i][2]

for i in range(len(s8)):
    s8[i][0] = s7[i][0]
    s8[i][1] = s7[i][1] * -1
    s8[i][2] = s7[i][2]

sphere1 = s7[s7[:, 1].argsort()]

number_of_particles1 = len(sphere1)

for n in range(number_of_particles1):
    #    if n < number_of_particles1 / 2 :
    X = np.asarray(sphere1[n], dtype=float)
    V = np.array([0, 0, 0], dtype=float)
    A = np.array([0, 0, 0], dtype=float)
    IntForce = np.array([0, 0, 0], dtype=float)
    ExtForce = np.array([0, 0, 0], dtype=float)
    radius = 0.1
    H = 0
    mass = 1
    density = 1
    volume = 1
    pressure = 1
    color = np.array([0, 0, 1, 1], dtype=float)
    RestPosition = np.array(sphere1[n], dtype=float)
    neighbors = []
    Plastic_epsilon = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]], dtype=float)
    particle = pt.Particle(env, X, V, A, IntForce, ExtForce, radius, H, mass, density, volume, pressure, color,
                           RestPosition, neighbors, Plastic_epsilon)
    # pf.Particle.addForce(particle)
    env.addParticle2Sphere1(particle)
    env.addParticle(particle)

# fixing some particles
# for i in range(6, number_of_particles1):
#     env.Sphere1[i].IsFixed()

# finding kernel radius
env.finding_Kernel_radius_Sphere1(sphere1, 15)

# finding nearest neighbours
env.NeighborsOfSphere1(sphere1)


# computing mass
env.Compute_ParticleMassOfSphere1(s=ph.s, rho_uni=ph.ro_uni)

# computing particle density
env.Compute_densityOfSphere1()

# computing particle volume
env.Compute_volumeOfSphere1()

# apply color map
env.ColorMapOfSphere1()

# compute inverse Moment matrix
for i, p1 in enumerate(env.Sphere1):
    env.Inv_MomentMatrix(p1)


def applyExtForce2Sphere1(env):
    for p1 in env.particles[0:int((len(env.particles) / 2)) + 1]:
        for p2 in env.particles[int((len(env.particles) / 2) + 1):len(env.particles)]:
            p1.attract(p2)

sphere2 = s8[s8[:, 1].argsort()]
number_of_particles2 = len(sphere2)
for n in range(number_of_particles2):
    #    if n < number_of_particles1 / 2 :
    X = np.asarray(sphere2[n], dtype=float)
    V = np.array([0, 0, 0], dtype=float)
    A = np.array([0, 0, 0], dtype=float)
    IntForce = np.array([0, 0, 0], dtype=float)
    ExtForce = np.array([0, 0, 0], dtype=float)
    radius = 0.1
    H = 0
    mass = 1
    density = 1
    volume = 1
    pressure = 1
    color = np.array([0, 0, 1, 1], dtype=float)
    RestPosition = np.array(sphere2[n], dtype=float)
    neighbors = []
    Plastic_epsilon = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]], dtype=float)
    particle = pt.Particle(env, X, V, A, IntForce, ExtForce, radius, H, mass, density, volume, pressure, color,
                           RestPosition, neighbors, Plastic_epsilon)
    # pf.Particle.addForce(particle)
    env.addParticle2Sphere2(particle)

# fixing some particles
# for i in range(1):
#     env.Sphere2[i].IsFixed()

# finding kernel radius
env.finding_Kernel_radius_Sphere2(sphere2, 15)

# finding nearest neighbours
env.NeighborsOfSphere2(sphere2)


# computing mass
env.Compute_ParticleMassOfSphere2(s=ph.s, rho_uni=ph.ro_uni)

# computing particle density
env.Compute_densityOfSphere2()

# computing particle volume
env.Compute_volumeOfSphere2()

# apply color map
env.ColorMapOfSphere2()

# compute inverse Moment matrix
for i, p1 in enumerate(env.Sphere2):
    env.Inv_MomentMatrix(p1)

env.particles = env.Sphere1 + env.Sphere2

def applyExtForce2Sphere2(env):
    for p1 in env.particles[int((len(env.particles)/2)+1):len(env.particles)]:
        for p2 in env.particles[0:int((len(env.particles)/2))+1]:
            p1.attract(p2)


