#!/usr/bin/env python3

import imageio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Import the 3D plotting toolkit

class Plotter:
    def __init__(self, simulation):
        self.simulation = simulation

    def plot_energy(self):
        energies = self.simulation.energies
        steps = range(len(energies))
        
        plt.figure(figsize=(10, 6))
        plt.plot(steps, energies, marker='.', linestyle='-', color='b', markersize=2)
        plt.title('Energy Variation Over Time')
        plt.xlabel('Step')
        plt.ylabel('Energy')
        plt.grid(True)
        plt.savefig('energy_over_time.png')
        plt.close()

    def plot_xz_projection(self, step_index=-1):
        system_snapshot = self.simulation.steps[step_index]
        xs, _, zs = system_snapshot
        max_x, _, max_z = self.simulation.container.corner

        plt.figure(figsize=(10, 6))
        plt.scatter(xs, zs, color='b')
        plt.xlim(0, max_x)
        plt.ylim(0, max_z)
        plt.title('Projection of the System on the XZ Plane')
        plt.xlabel('X')
        plt.ylabel('Z')
        plt.grid(True)
        plt.savefig(f'xz_projection_step_{step_index}.png')
        plt.close()

    def plot_3d_state(self, step_index=-1):
        """
        Plots the 3D state of the system for a given snapshot.
        If step_index is not provided, it plots the last snapshot by default.
        """
        system_snapshot = self.simulation.steps[step_index]
        xs, ys, zs = system_snapshot  # Extract xs, ys, and zs from the snapshot

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')  # Create a 3D subplot
        ax.scatter(xs, ys, zs, color='g', marker='o')  # Plot the points in 3D space

        ax.set_title('3D State of the System')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.grid(True)
        plt.savefig(f'3d_state_step_{step_index}.png')
        plt.close()
    
    def create_gif(self, filename='system_evolution.gif', plot_type='3d', step_interval=100):
        """
        Creates a GIF showing the evolution of the system over time in either 2D or 3D.

        Parameters:
        - filename: Name of the output GIF file.
        - plot_type: Type of plot ('2d' or '3d') for the GIF.
        - step_interval: Interval between steps to include in the GIF.
        """
        images = []
        for step_index in range(0, len(self.simulation.steps), step_interval):
            fig = plt.figure(figsize=(10, 8))
            if plot_type == '3d':
                ax = fig.add_subplot(111, projection='3d')
                system_snapshot = self.simulation.steps[step_index]
                xs, ys, zs = system_snapshot
                ax.scatter(xs, ys, zs, color='b', marker='o')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
            elif plot_type == '2d':
                ax = fig.add_subplot(111)
                system_snapshot = self.simulation.steps[step_index]
                xs, _, zs = system_snapshot  # Assuming XZ projection for 2D
                ax.scatter(xs, zs, color='b')
                ax.set_xlabel('X')
                ax.set_ylabel('Z')
            else:
                raise ValueError("plot_type must be '2d' or '3d'")
            
            ax.set_title(f'Step {step_index}')
            # Save plot to a temporary PNG file
            temp_filename = f'temp_step_{step_index}.png'
            plt.savefig(temp_filename)
            plt.close()
            # Append the image to the list of images
            images.append(imageio.imread(temp_filename))
            # Optionally, remove the temporary file here if desired
        
        # Create GIF
        imageio.mimsave(filename, images, fps=10)