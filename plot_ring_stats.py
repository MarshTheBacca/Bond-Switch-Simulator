import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

cwd = Path(__file__).parent

with open(cwd.joinpath("testing/output_files/test_ringstats.out"), "r") as file:
    ring_stats = [[i, [float(n) for n in line.rstrip().split()]] for i, line in enumerate(file.readlines())]
ring_stats = np.array([n for _, n in ring_stats])
for i in range(ring_stats.shape[1]):
    if np.count_nonzero(ring_stats[:,i]) > 0:
        plt.plot(ring_stats[:,i], label="Ring size " + str(i+3))

plt.title("Pore, 100 thermalisation steps, 4000 annealing steps, switch acceptance 11.7%")
plt.xlabel("Step")
plt.ylabel("Proportion of rings of size n")
plt.legend()
plt.show()