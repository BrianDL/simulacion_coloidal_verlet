#!/usr/bin/env python3

import numpy as np
from plotter import Plotter
from container import Region, get_region

# Modify the Simulation class to include energies
class Simulation:
    def __init__(self, steps, L, energies=None):
        self.steps = steps
        self.container = get_region('block', (L,L,L))
        self.energies = energies

def read_energias(filename):
    energies = []
    with open(filename, 'r') as f:
        for line in f:
            energies.append(float(line.strip()))
    return np.array(energies)

# Modify the read_estados function to also read energies
def read_estados(estados_filename, energias_filename, L=1000):
    steps = []
    with open(estados_filename, 'r') as f:
        while True:
            x = f.readline().strip()

            if not x: break

            y = f.readline().strip()
            z = f.readline().strip()

            step = (
                np.array(list(map(float, x.split()))),
                np.array(list(map(float, y.split()))),
                np.array(list(map(float, z.split())))
            )
            steps.append(step)

    energies = read_energias(energias_filename)

    return Simulation(steps, L, energies)
# Read the data
sim = read_estados('estados.txt', 'energias.txt')

# You can now access the energies like this:
print(f"Number of energy values: {len(sim.energies)}")
print(f"First energy value: {sim.energies[0]}")
print(f"Last energy value: {sim.energies[-1]}")