import shutil
import numpy as np
from pathlib import Path
from typing import Callable, Tuple


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


def get_nums(path: Path, num_atoms_line=3, num_bonds_line=5) -> Tuple[int, int]:
    num_atoms = int(get_line(path, num_atoms_line).split()[0])
    num_bonds = int(get_line(path, num_bonds_line).split()[0])
    return num_atoms, num_bonds


def write_net(path: Path, bonds: np.array, num_atoms: int):
    # It may appear {prefix}_net.dat is erroneous for Si2O3 as it does not seem to be atom 0 at the top of the file,
    # but I have triple checked this is not an error, and simply due to the fact the results.dat file does not have
    # atoms sorted by ID, and the list of bonds does not identify all bonds to one node and move on to the next node,
    # as I can see nodes already listed way further down in the list bonded to other nodes.

    nodes = [[] for atom in np.arange(num_atoms)]
    for bond in bonds:
        selected_node = bond[2] - 1  # Because NetMC IDs atoms from 0 and LAMMPS doesnt
        bonded_node = bond[3] - 1  # Because NetMC IDs atoms from 0 and LAMMPS doesnt
        nodes[selected_node].append(bonded_node)
        nodes[bonded_node].append(selected_node)
    #           atom 0          atom 1 ...
    # nodes = [ [n1, n2, n3], [n4, n5, n6], ...]
    with open(path, "w") as net_file:
        for node in nodes:
            line = ""
            for bonded_node in node:
                line += f"{bonded_node:<10}"
            net_file.write(f"{line}\n")


def write_coords(path: Path, coords: np.array):
    with open(path, "w") as coords_file:
        for coord in coords:
            coords_file.write(f"{coord[0]:<24}{coord[1]:<24}{coord[2]:<24}\n")


def write_aux(path: Path, num_atoms, dims):
    with open(path, "w") as aux_file:
        aux_file.write(f"{num_atoms}\n")
        aux_file.write(f"{3:<10}{3:<10}\n")
        aux_file.write("2DE\n")
        aux_file.write(f"{dims[0]:<24.6f}{dims[1]:<24.6f}\n")
        aux_file.write(f"{1 / dims[0]:<24.6f}{1 / dims[1]:<24.6f}\n")


def process_results(results_folder: Path, results_file_name: str, netmc_prefix: str):
    results_file_path = results_folder.joinpath(results_file_name)
    num_atoms, num_bonds = get_nums(results_file_path)

    # Identify the line at which dimensions are listed, ie, where we see 'xlo xhi'
    skip_to_dims = (get_target_line(results_file_path, lambda line: line.endswith("xlo xhi")) - 1)

    # Skip to this line and read the next 3 into an array, taking the second column for maximums for our axes
    # z axis maximum is not actually used in this program as it is not written to the aux file
    dims = np.genfromtxt(results_file_path, skip_header=skip_to_dims, max_rows=3)[:, 1]

    # Identify the line at which atoms are listed, ie, the line after the file says "Atoms"
    skip_to_atoms = (get_target_line(results_file_path, lambda line: line.startswith("Atoms")) + 1)

    # Generate coordinates array for atoms, x -> col4, y -> col5, z -> col6, then export to new NetMC file
    atoms_array = np.genfromtxt(results_file_path, skip_header=skip_to_atoms, max_rows=num_atoms)
    coords = atoms_array[:, 3:6].astype(np.float64)
    write_coords(results_folder.joinpath(f"{netmc_prefix}_crds.dat"), coords)

    # There are two titles (3 lines each) and 2 sets of num_atoms lines (coords and velocities) before bonds
    # This can be a dangerous assumption for finding where bonds are listed!!
    bonds = np.genfromtxt(results_file_path,
                          skip_header=skip_to_atoms + 2 * num_atoms + 6,
                          max_rows=num_bonds,
                          dtype=np.int64)
    write_net(results_folder.joinpath(f"{netmc_prefix}_net.dat"), bonds, num_atoms)

    # Copy over dual.dat and write the aux file
    shutil.copyfile(results_folder.joinpath("test_A_dual.dat"),
                    results_folder.joinpath(f"{netmc_prefix}_dual.dat"))
    write_aux(results_folder.joinpath(f"{netmc_prefix}_aux.dat"), num_atoms, dims)


cwd = Path(__file__).parent

# Find initial and final folders from netmc.inpt
netmc_input_path = cwd.joinpath("netmc.inpt")

# initial_folder = cwd.joinpath(Path(get_line(netmc_input_path, 6).split()[0][2:])) #(Currently 'Results', unused)
final_folder = cwd.joinpath(Path(get_line(netmc_input_path, 2).split()[0][2:]))  # (Currently 'Run')

si_results_path = final_folder.joinpath("Si_results.dat")
if si_results_path.is_file():
    print("Writing Si Positions\n")
    process_results(final_folder, "Si_results.dat", "test_Si")

bn_results_path = final_folder.joinpath("BN_results.dat")
if bn_results_path.is_file():
    print("Writing BN Positions\n")
    process_results(final_folder, "BN_results.dat", "test_BN")

si2o3_results_path = final_folder.joinpath("Si2O3_results.dat")
if si2o3_results_path.is_file():
    print("Writing Si2O3 Positions\n")
    process_results(final_folder, "Si2O3_results.dat", "test_Si2O3")
