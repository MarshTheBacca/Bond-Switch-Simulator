from __future__ import annotations
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np



@dataclass
class NetMCOutputData:
    file_path: Path

    def __post_init__(self):
        self.steps = []
        self.temperatures = [] 
        self.energies = []
        self.entropies = []
        self.pearson_coeffs = []
        self.aboav_weaires = []
        self.ring_sizes = []
        self.ring_areas = []
        with open(self.file_path, "r") as file:
            for _ in range(4):
                file.readline()
            lines = file.readlines()
        for line in lines:
            line = line.split(",")
            self.steps.append(int(line[0]))
            self.temperatures.append(float(line[1]))
            self.energies.append(float(line[2]))
            self.entropies.append(float(line[3]))
            self.pearson_coeffs.append(float(line[4]))
            self.aboav_weaires.append(float(line[5]))
            self.ring_sizes.append({int(ring_size): float(proportion) for ring_size, proportion in (pair.split(":") for pair in line[6].split(";"))})
            self.ring_areas.append([float(area) for area in line[7].split(";")])
            
    def plot_energy(self):
        plt.plot(self.steps, self.energies)
        plt.title("Energy change over time")
        plt.xlabel("Step")
        plt.ylabel("Energy (Hartrees)")
        plt.show()
    
    def plot_temperature(self):
        plt.plot(self.steps, self.temperatures)
        plt.title("Temperature change over time")
        plt.xlabel("Step")
        plt.ylabel("Temperature (au)")
        plt.show()
    
    def plot_entropy(self):
        plt.plot(self.steps, self.entropies)
        plt.title("Entropy change over time")
        plt.xlabel("Step")
        plt.ylabel("Entropy (Hartrees Temperature^-1)")
        plt.show()
    
    def plot_pearsons_coeffs(self):
        plt.plot(self.steps, self.pearson_coeffs)
        plt.ylim(-1, 1)
        plt.title("Pearson's Correlation Coefficient change over time")
        plt.xlabel("Step")
        plt.ylabel("Pearson's Correlation Coefficient")
        plt.show()
    
    def plot_aboav_weaire(self):
        plt.plot(self.steps, self.aboav_weaires)
        plt.title("Aboav-Weaire parameter change over time")
        plt.xlabel("Step")
        plt.ylabel("Aboav-Weaire parameter")
        plt.show()
        
    def plot_ring_sizes(self):
        # Get all unique ring sizes across all steps
        all_ring_sizes = sorted(set(k for d in self.ring_sizes for k in d.keys()))

        # Create a 2D list where each sublist represents a ring size over time
        ring_sizes_transposed = [[step.get(ring_size, 0) for step in self.ring_sizes] for ring_size in all_ring_sizes]

        # Filter out ring sizes that are always zero
        ring_sizes_filtered = [ring for ring in ring_sizes_transposed if sum(ring) > 0]

        # Create labels for the ring sizes
        ring_labels = ['Ring Size ' + str(ring_size) for ring_size, ring in zip(all_ring_sizes, ring_sizes_transposed) if sum(ring) > 0]

        plt.stackplot(self.steps, *ring_sizes_filtered, labels=ring_labels)
        plt.title('Proportion of Ring Sizes per Step')
        plt.xlabel('Step')
        plt.ylabel('Proportion')

        # Set the limits of the x-axis and y-axis
        plt.gca().set_xlim([0, max(self.steps)])
        plt.gca().set_ylim([0, 1])

        # Place the legend outside of the plot
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.show()
            
    def plot_ring_areas(self):
        plt.plot(self.steps, [np.mean(areas) for areas in self.ring_areas])
        plt.title("Average Ring Area change over time")
        plt.xlabel("Step")
        plt.ylabel("Average Area (Bohr Radii^2)")
        plt.show()
        
        
            
    