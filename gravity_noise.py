import numpy as np
import math
import csv
import random
import pandas as pd

###
total_time = 1e5
spawn_time = 1e4
dt = 0.0001
###

Mpc = 5.38552341e20
seconds_year = 365.25 * 24 * 3600
G = 6.6743e-11 * seconds_year ** 2
c = 3e8 * seconds_year
r = 100 * Mpc
Msun = 1.989e30
num = int(total_time / dt)
T = 24 * 3600 
gsi = 6.6743e-11
total_mass = 2 * 1.4 * Msun
dist0 = ((T ** 2 * gsi * total_mass) / (2 * np.pi ** 2))**(1/3)
filename = 'output.csv'
chunks = int(spawn_time/dt)

class Particle():
    def __init__(self, x, y, density, radius, color, omega):
        self.x = x
        self.y = y
        self.density = density
        self.radius = radius
        self.mass = density * (4/3) * np.pi * radius ** 3
        self.color = color
        self.omega = omega

        self.x_vel = 0
        self.y_vel = 0

def gravity(particles, G, dt):
    if len(particles) == 2:
        p1 = particles[0]
        p2 = particles[1]

        dx = p1.x - p2.x
        dy = p1.y - p2.y

        dist2 = dx**2 + dy**2
        dist = np.sqrt(dist2)
        if dist < p1.radius + p2.radius:
            return True
        F = (G * p1.mass * p2.mass) / dist2
        Fx = F * dx / dist
        Fy = F * dy / dist

        p1.x_vel -= Fx / p1.mass * dt
        p1.y_vel -= Fy / p1.mass * dt

        p2.x_vel += Fx / p2.mass * dt
        p2.y_vel += Fy / p2.mass * dt

        LGW = (32 * G ** 4 * (p1.mass + p2.mass) ** 3 * ((p1.mass * p2.mass)/(p1.mass + p2.mass))** 2)/(5 * c ** 5 * (dist/2) ** 5)
        
        Eloss = LGW * dt
        Ecurr = -G * p1.mass * p2.mass / (dist)
        ENew = Ecurr - Eloss
        rnew = -G * p1.mass * p2.mass / (2 * ENew)
        scale = rnew /(dist/2)     
        #print(scale)
        p1.x = masscentrx + (p1.x - masscentrx) * scale
        p1.y = masscentry + (p1.y - masscentry) * scale

        p2.x = masscentrx + (p2.x - masscentrx) * scale
        p2.y = masscentry + (p2.y - masscentry) * scale

        scalev = 1 / np.sqrt(scale)
        p1.x_vel *= scalev
        p1.y_vel *= scalev
        p2.x_vel *= scalev
        p2.y_vel *= scalev

        p1.omega = np.sqrt(G * (p1.mass + p2.mass) / (dist/2)**3)
        p2.omega = p1.omega

    return False

masscentrx = 0
masscentry = 0

def GW(omega, t, R, r, m):
    Amp = (8 * m * omega ** 2 * R ** 2)/r 
    hplusz = -Amp * np.cos(2 * omega * t)
    hplusx = hplusz/2
    hxz = Amp * np.sin(2 * omega * t)
    hxx = hxz/2
    return hplusz, hxz, hplusx, hxx

def datacol(my_particles, data1, data2, data3, data4):
    t = 0
    j = 0

    running = True
    while running and j < chunks:
        if gravity(my_particles, G, dt):
            running = False
        
        hpz, hxz, hpx, hxx = GW(star1.omega, t, (star1.x - masscentrx), r, star1.mass)
        data1[j] += hpz
        data2[j] += hxz
        data3[j] += hpx
        data4[j] += hxx

        for p in my_particles:
            p.x += p.x_vel * dt
            p.y += p.y_vel * dt
        
        t += dt
        j += 1
    return data1, data2, data3, data4

iter = math.floor(total_time/spawn_time)
print("Setup done")

total_particles = []
for i in range(iter):
    theta = random.uniform(0, 2 * np.pi)
    star1 = Particle(-np.sin(theta) * (dist0/2), np.cos(theta) * (dist0/2), (1.4 * Msun * 3)/(4 * np.pi * 14000 ** 3), 14000, (255, 255, 255), 0)
    star2 = Particle(-star1.x, -star1.y, (1.4 * Msun * 3)/(4 * np.pi * 14000 ** 3), 14000, (255, 255, 255), 0)
    theta = random.uniform(0, 2 * np.pi)
    v0 = np.sqrt((G * star1.mass) / (4 * (dist0)/2))
    star1.x_vel = -v0 * np.cos(theta)
    star1.y_vel = v0 * np.sin(theta)
    star2.x_vel = -star1.x_vel
    star2.y_vel = -star1.y_vel
    my_particles = [star1, star2]
    total_particles.append(my_particles)

with open(filename, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    for i in range(int(total_time/spawn_time)):
        print(f"Starting chunk {i+1}")
        z = 1
        data1 = [0.0] * chunks
        data2 = [0.0] * chunks
        data3 = [0.0] * chunks
        data4 = [0.0] * chunks
        for j in range(0, z):
            print(f"Starting merger {j+1}")
            data1, data2, data3, data4 = datacol(total_particles[j], data1, data2, data3, data4)
            df = zip(data1, data2, data3, data4)
            print("Outputting data")
            writer.writerows(df)
        z += 1

print("All done")