#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def read_estados(filename):
    states = []
    with open(filename, 'r') as f:
        while True:
            x = f.readline().strip()
            if not x:
                break
            y = f.readline().strip()
            z = f.readline().strip()
            
            state = np.array([
                list(map(float, x.split())),
                list(map(float, y.split())),
                list(map(float, z.split()))
            ])
            states.append(state)
    return states

# Read the data
states = read_estados('estados.txt')

# Set up the figure and 3D axis
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Initialize the scatter plot
scat = ax.scatter([], [], [])

# Animation update function
def update(frame):
    ax.clear()
    state = states[frame]
    scat = ax.scatter(state[0], state[1], state[2], c='b', s=50)
    ax.set_xlim(min(state[0]), max(state[0]))
    ax.set_ylim(min(state[1]), max(state[1]))
    ax.set_zlim(min(state[2]), max(state[2]))
    ax.set_title(f'Frame {frame}')
    return scat,

# Create the animation
anim = FuncAnimation(fig, update, frames=len(states), interval=200, blit=False)

# Display the animation
# plt.show()

# Optionally, save the animation (uncomment to use)
anim.save('system_evolution.gif', writer='pillow', fps=5)