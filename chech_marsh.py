import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import colormaps as cm
from pathlib import Path
from typing import Callable

# Question for Olly: if we read from a Si2O3 coordinates file, are the first 2/5 of coordinates
# for Si atoms, and the next 3/5 oxygen atoms?
# What about SiO2?


def get_graph(network: dict) -> nx.Graph:
    network_graph = nx.Graph()
    for selected_node in network:
        network_graph.add_node(selected_node, pos=network[selected_node]["crds"])
        for bonded_node in network[selected_node]["net"]:
            network_graph.add_edge(selected_node, bonded_node)
    return network_graph


def ordered_cnxs(network_A: nx.Graph, network_B: nx.Graph):
    for ring in network_B:
        cnxs = network_B[ring]["dual"]
        new_cnxs, new_coords = [], []
        new_cnxs.append(cnxs[0])
        new_coords.append(network_A[cnxs[0]]["crds"])
        i = 0
        while i < len(cnxs) - 1:
            node0 = new_cnxs[i]
            connected_to_0 = network_A[node0]["net"]
            options = set(connected_to_0).intersection(cnxs)
            for val in new_cnxs:
                if val in options:
                    options.remove(val)
            options = list(options)
            new_cnxs.append(options[0])
            new_coords.append(network_A[options[0]]["crds"])
            i += 1
        area = 0
        for i in range(len(cnxs)):
            x0, y0, x1, y1 = (new_coords[i - 1][0],
                              new_coords[i - 1][1],
                              new_coords[i][0],
                              new_coords[i][1])
            area += x0 * y1 - x1 * y0
        if area > 0:
            new_cnxs.reverse()
        network_B[ring]["dual"] = new_cnxs


def get_line(path: Path, target_line: int) -> str:
    with open(path, "r") as file:
        for line_num, line in enumerate(file):
            if line_num + 1 == target_line:
                return line


def get_target_line(path: Path, func: Callable[[str], bool]) -> int:
    with open(path, "r") as file:
        for line_num, line in enumerate(file):
            if func(line.rstrip()):
                return line_num + 1


def read_aux(path: Path):
    num_atoms = int(get_line(path, 1))
    xlo = float(get_line(path, 4).split()[0])
    xhi = float(get_line(path, 4).split()[1])
    ylo = float(get_line(path, 5).split()[0])
    yhi = float(get_line(path, 5).split()[1])
    return num_atoms, xlo, xhi, ylo, yhi


def read_crds_net_dual_aux(path: Path, prefix: str, file_type, prefixes=None):
    if not prefixes:
        prefixes = [prefix] * 4
    aux_data = read_aux(path.joinpath(f"{prefixes[0]}_aux.dat"))
    with open(path.joinpath(f"{prefixes[1]}_crds.dat"), "r") as crds_file:
        coords = np.genfromtxt(crds_file, dtype=np.float64)
    if file_type == "Si2O3":
        num_nodes = int(aux_data[0] * 2 / 5)
        oxygen_coords = coords[num_nodes:]
        coords = coords[:num_nodes]
    elif file_type == "SiO2":
        num_nodes = int(aux_data[0] / 6)
        oxygen_coords = coords[num_nodes:]
        coords = coords[:num_nodes]
    else:
        num_nodes = aux_data[0]
        oxygen_coords = None
    net, dual = [], []  # Because nodes file is not homogenous
    with open(path.joinpath(f"{prefixes[2]}_net.dat"), "r") as net_file:
        for line in net_file:
            net.append(np.array(line.strip().split(), dtype=np.int64))
    with open(path.joinpath(f"{prefixes[3]}_dual.dat"), "r") as dual_file:
        for line in dual_file:
            dual.append(np.array(line.strip().split(), dtype=np.int64))

    return coords, net, dual, num_nodes, aux_data[1], aux_data[2], oxygen_coords


ring_colours = ["white", "white", "white", cm["GnBu"](140), cm["Greens"](100), cm["Blues"](150),
                cm["Greys"](90), cm["Reds"](105), cm["YlOrBr"](100), cm["PuRd"](100),
                cm["RdPu"](80)] + 40 * ["black"]

cwd = Path(__file__).parent
folder_name = "Run"
file_name = "test"
file_type = "Si2O3"
folder_path = cwd.joinpath(folder_name)

coords_A, net_A, dual_A, num_nodes_A, xlo_A, xhi_A, oxygen_coords_A = read_crds_net_dual_aux(folder_path, f"{file_name}_A", "unknown",
                                                                                             # prefixes=["test_Si2O3", "test_Si2O3", "test_A", "test_A"]
                                                                                             )
coords_B, net_B, dual_B, num_nodes_B, xlo_B, xhi_B, oxygen_coords_A = read_crds_net_dual_aux(folder_path, f"{file_name}_B", "unknown"
                                                                                             # prefixes=["test_A", "test_B", "test_B", "test_B"]
                                                                                             )
coords_B = np.multiply(coords_B, xlo_A / xlo_A)

print("*" * 40, f"Number of network_A nodes: {num_nodes_A}", "*" * 40)
print("*" * 40, f"Number of network_B nodes: {num_nodes_B}", "*" * 40)

network_A = {i: {"crds": coords_A[i], "net": net_A[i], "dual": dual_A[i]} for i in np.arange(num_nodes_A)}
network_B = {i: {"crds": coords_B[i], "net": net_B[i], "dual": dual_B[i]} for i in np.arange(num_nodes_B)}

# network_X = {0:{"crds": [coords_0], "net": [net_0], "dual": [dual_0]}, 1: {... }}

ordered_cnxs(network_B, network_A)
broken_nodes = []

print("Checking for broken nodes...\n")
for selected_node in network_A:
    if selected_node not in network_A:
        print(f"Node not found in network_A: {selected_node}")
    for bonded_node in network_A[selected_node]["net"]:
        if bonded_node not in network_A:
            print(f"Bonded node not found in network_A: {bonded_node}")
        if selected_node not in network_A[bonded_node]["net"]:
            print("#" * 40)
            print("Selected node not found in bonded node's network")
            print(f"selected_node: {selected_node}\tbonded_node: {bonded_node}")
            print(f"selected_node['net']: {selected_node['net']}\tbonded_node['net']: {bonded_node['net']}\n")
            broken_nodes.append(selected_node)
            broken_nodes.append(bonded_node)
if len(broken_nodes) == 0:
    print("No broken nodes detected!\n")

for node in network_A:
    for ring in network_A[node]["dual"]:
        if node not in network_B[ring]["dual"]:
            print("Broken Node-Dual")

network_A_graph, network_B_graph = get_graph(network_A), get_graph(network_B)

print("Plotting 1")
plt.scatter(coords_A[:, 0], coords_A[:, 1], color="b")
plt.scatter(coords_B[:, 0], coords_B[:, 1], color="g")
for node in broken_nodes:
    plt.text(coords_A[node, 0], coords_A[node, 1], node, fontsize=8)
for selected_node in network_A:
    for bonded_node in network_A[selected_node]["net"]:
        dx = abs(network_A[selected_node]["crds"][0] - network_A[bonded_node]["crds"][0])
        dy = abs(network_A[selected_node]["crds"][1] - network_A[bonded_node]["crds"][1])
        if dx < 20 and dy < 20:
            plt.plot([network_A[selected_node]["crds"][0], network_A[bonded_node]["crds"][0]],
                     [network_A[selected_node]["crds"][1], network_A[bonded_node]["crds"][1]],
                     color="b", alpha=0.5)
plt.show()

print("Plotting 2")
fig = plt.figure()
ax = fig.add_subplot(111)
patches, plotting_ring_colours = [], []

for selected_node in network_B:
    ring_atoms = network_B[selected_node]["dual"]
    ring_atom_coords = [network_A[ring_atoms[0]]["crds"][:2]]
    for i in range(1, len(ring_atoms)):
        dx = network_A[ring_atoms[i]]["crds"][0] - ring_atom_coords[i - 1][0]
        dy = network_A[ring_atoms[i]]["crds"][1] - ring_atom_coords[i - 1][1]
        if dx > xlo_A / 2:
            dx -= xlo_A
        elif dx < -xlo_A / 2:
            dx += xlo_A
        if dy > xhi_A / 2:
            dy -= xhi_A
        elif dy < -xhi_A / 2:
            dy += xhi_A
        ring_atom_coords.append([network_A[ring_atoms[i]]["crds"][0] + dx, ring_atom_coords[i - 1][0] + dy])
    ring_atom_coords = np.array(ring_atom_coords)
    for x in range(0, 2):
        for y in range(0, 2):
            new_ring_atom_coords = ring_atom_coords.copy()
            for i in range(ring_atom_coords.shape[0]):
                new_ring_atom_coords[i, 0] += x * xlo_A
                new_ring_atom_coords[i, 1] += y * xhi_A
            patches.append(Polygon(new_ring_atom_coords, closed=True))
            plotting_ring_colours.append(ring_colours[len(new_ring_atom_coords)])

ax.add_collection(PatchCollection(patches,
                                  facecolor=plotting_ring_colours,
                                  edgecolor="k",
                                  linewidths=0.5,
                                  alpha=0.5,
                                  zorder=1))

if oxygen_coords_A:
    plt.scatter(oxygen_coords_A[:, 0], oxygen_coords_A[:, 1], color="r", s=0.5)

ax.set_xlim(-xlo_A * 0.5, xlo_A * 2.5)
ax.set_ylim(-xhi_A * 0.5, xhi_A * 2.5)
plt.show()

print("Plotting 3")

for selected_node in network_A:
    selected_node_dual = network_A[selected_node]["dual"]
    for i in range(len(selected_node_dual)):
        dx = abs(network_B[selected_node_dual[i - 1]]["crds"][0] - network_B[selected_node_dual[i]]["crds"][0])
        dy = abs(network_B[selected_node_dual[i - 1]]["crds"][1] - network_B[selected_node_dual[i]]["crds"][1])
        if dx < 20 and dy < 20:
            plt.plot([network_B[selected_node_dual[i - 1]]["crds"][0], network_B[selected_node_dual[i]]["crds"][0]],
                     [network_B[selected_node_dual[i - 1]]["crds"][1], network_B[selected_node_dual[i]]["crds"][1]],
                     alpha=0.5)
plt.show()
