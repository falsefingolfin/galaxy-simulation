import pandas as pd
import numpy as np


# NEED TO ADD TILT TO GALAXIES


rmin = 25.0                         # kpc
G = 1                               # (kpc)^3 / ((M_sol)(yr * 10^8)^2)
M = 10**11 * (4.498 * 10**(-6))     # 10^11 solar masses in G=1 units)
m = M * 10**(-12)
ecc = 0.5
wa = np.radians(-30)
ia = np.radians(60)
wb = np.radians(30)
ib = np.radians(-60)


# tilt rotation matrices
rotxa = np.array([[1, 0, 0],
                 [0, np.cos(wa), -1*np.sin(wa)],
                 [0, np.sin(wa), np.cos(wa)]])

rotya = np.array([[np.cos(ia), 0, np.sin(ia)],
                 [0, 1, 0],
                 [-1*np.sin(ia), 0, np.cos(ia)]])

rotza = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])   
temprot = np.matmul(rotza, rotya)
rotmatrixa = np.matmul(temprot, rotxa)


rotxb = np.array([[1, 0, 0],
                 [0, np.cos(wb), -1*np.sin(wb)],
                 [0, np.sin(wb), np.cos(wb)]])

rotyb = np.array([[np.cos(ib), 0, np.sin(ib)],
                 [0, 1, 0],
                 [-1*np.sin(ib), 0, np.cos(ib)]])

rotzb = np.array([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]])       
temprot = np.matmul(rotzb, rotyb)
rotmatrixb = np.matmul(temprot, rotxb)


# create galaxy of 345 disk particles, each particle has 6 parameters [x,y,z,vx,vy,vz,m]
# leftGalaxy[0] is the center of the galaxy point, everything else is a disk particle
leftGalaxyPos = np.zeros(shape=(342,3))
leftGalaxyVel = np.zeros(shape=(342,3))
rightGalaxyPos = np.zeros(shape=(342,3))
rightGalaxyVel = np.zeros(shape=(342,3))

leftGalaxy = np.zeros(6)
rightGalaxy = np.zeros(6)


def populateGalaxy(positions, velocities, orientation, rotmatrix):
    # rotation = 1: counterclockwise        rotation = -1, clockwise
    e = 0.2*rmin     # softening
    n = 12
    r = 0.2*rmin
    index = 0

    for k in range(12):
        v0 = np.sqrt(G*M*r / (r**2 + e**2))       # v0 for circular orbit

        # fill each ring with n particles
        for j in range(n):
            angle = j * 2*np.pi / n
            positions[index][0] = r * np.cos(angle)
            positions[index][1] = r * np.sin(angle)
            positions[index][2] = 0
            velocities[index][0] = orientation * v0 * np.cos(angle + np.pi/2)
            velocities[index][1] = orientation * v0 * np.sin(angle + np.pi/2)
            velocities[index][2] = 0

            # rotate all positions and velocities by rotation matrix
            positions[index] = rotmatrix @ positions[index]
            velocities[index] = rotmatrix @ velocities[index]

            index += 1

        n += 3
        r += 0.05*rmin


def shiftCenter(galaxy, positions, velocities, direction):
    # direction = 1: left galaxy        direction = -1: right galaxy
    a = rmin/2 / (1-ecc)
    # r = rmin * (1+e) / (1 - e * cos(0)) = a(1+e)
    r = a * (1+ecc)
    v = np.sqrt(G*M * (2/r - 1/a)) / 2
    galaxy[0] = -1 * direction * r
    galaxy[4] = direction * v

    positions[:,0] += galaxy[0]
    velocities[:,1] += galaxy[4]



populateGalaxy(leftGalaxyPos, leftGalaxyVel, -1, rotmatrixb)
populateGalaxy(rightGalaxyPos, rightGalaxyVel, -1, rotmatrixa)
shiftCenter(leftGalaxy, leftGalaxyPos, leftGalaxyVel, 1)
shiftCenter(rightGalaxy, rightGalaxyPos, rightGalaxyVel, -1)



aarsethParams = "{n} {eta} {dt} {T} {e}\n"         # {number of particles} {integration timestep} {major output timestep} {total time} {smoothing}
particle = "{mass} {x:0.6f} {y:0.6f} {z:0.6f} {vx:0.6f} {vy:0.6f} {vz:0.6f}\n"

f = open("initc.data", "w")

# setup aarseth code parameters
f.write(aarsethParams.format(n=686, eta=0.01, dt=0.01, T=10.0, e=(0.2*rmin)**2))

# write galaxy centers to file
f.write(particle.format(mass=M, x=leftGalaxy[0], y=leftGalaxy[1], z=leftGalaxy[2],
                            vx=leftGalaxy[3], vy=leftGalaxy[4], vz=leftGalaxy[5]))
f.write(particle.format(mass=M, x=rightGalaxy[0], y=rightGalaxy[1], z=rightGalaxy[2],
                            vx=rightGalaxy[3], vy=rightGalaxy[4], vz=rightGalaxy[5]))

# write galaxy particles to file
for j in range(342):
    f.write(particle.format(mass=m, x=leftGalaxyPos[j][0], y=leftGalaxyPos[j][1], z=leftGalaxyPos[j][2],
                                vx=leftGalaxyVel[j][0], vy=leftGalaxyVel[j][1], vz=leftGalaxyVel[j][2]))
    f.write(particle.format(mass=m, x=rightGalaxyPos[j][0], y=rightGalaxyPos[j][1], z=rightGalaxyPos[j][2],
                                vx=rightGalaxyVel[j][0], vy=rightGalaxyVel[j][1], vz=rightGalaxyVel[j][2]))

f.close()