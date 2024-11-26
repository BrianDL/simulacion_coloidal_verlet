#!/usr/bin/env python3

import numpy as np
from plotter import Plotter
from container import Region, get_region

class Simulation:
    def __init__(self, steps, L):
        self.steps = steps
        self.container = get_region('block', (L,L,L) )

def read_estados(filename, L = 1000):
    steps = []
    with open(filename, 'r') as f:
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
    return Simulation(steps, L)

# Read the data
sim = read_estados('estados.txt')