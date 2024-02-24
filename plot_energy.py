import matplotlib.pyplot as plt
from pathlib import Path

cwd = Path(__file__).parent
energy = []
with open(cwd.joinpath("testing/output_files/test_energy.out"), "r") as file:
    data = [[i, float(line.rstrip())] for i, line in enumerate(file.readlines())]
    
plt.plot([i for i, _ in data], [energy for _, energy in data])
plt.title("Pore, 100 thermalisation steps, 4000 annealing steps, switch acceptance 11.7%")
plt.xlabel("Step")
plt.ylabel("Energy (Hartrees)")
plt.show()