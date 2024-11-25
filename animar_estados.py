#!/usr/bin/env python3

import numpy as np
from plotter import Plotter

class Simulation:
    def __init__(self, steps):
        self.steps = steps

def read_estados(filename):
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
    return Simulation(steps)

# Read the data
simulation = read_estados('estados.txt')