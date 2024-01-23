import shutil
from pathlib import Path


def get_line(path: Path, target_line: int) -> str:
    with open(path, "r") as file:
        for line_num, line in enumerate(file):
            if line_num + 1 == target_line:
                return line


def export_1D(path: Path, lines: list):
    string = ""
    for line in lines:
        string += line + "\n"
    with open(path, "w") as file:
        file.write(string)


def edit_line(path: Path, line: int, replacement) -> list:
    with open(path, "r") as file:
        lines = [line.rstrip("\n") for line in file]
    lines[line - 1] = replacement
    return lines


cwd = Path(__file__).parent

# Find initial and final folders from netmc.inpt
netmc_input_path = cwd.joinpath("netmc.inpt")

initial_folder = cwd.joinpath(Path(get_line(netmc_input_path, 6).split()[0][2:]))  # Results
final_folder = cwd.joinpath(Path(get_line(netmc_input_path, 2).split()[0][2:]))    # Run

# If the fixed_rings.dat exists in the inital folder, copy it over to the final folder.
# Else, create a new fixed_rings.dat in final folder with '0' inside
if initial_folder.joinpath("fixed_rings.dat").exists():
    shutil.copyfile(initial_folder.joinpath("fixed_rings.dat"), final_folder.joinpath("fixed_rings.dat"))
else:
    with open(final_folder.joinpath("fixed_rings.dat"), "w+") as fixed_rings_file:
        fixed_rings_file.write("0")

# Copy over the following files                                      These ones were commented out
files = ("PARM_C.lammps", "PARM_Si2O3.lammps", "PARM_Si.lammps")  # ("PARM_BN.lammps", "COULOMB.table")
for file in files:
    shutil.copyfile(initial_folder.joinpath(file), final_folder.joinpath(file))

# Copy over Si.in and replace line 15 with read_restart     Si_restart.restart
export_1D(path=final_folder.joinpath("Si.in"),
          lines=edit_line(path=initial_folder.joinpath("Si.in"),
                          line=15,
                          replacement=f"read_restart\t\t\t{Path(initial_folder.name).joinpath('Si_restart.restart')}"))

# Copy over Si.in and replace line 15 with read_restart     Si_restart.restart
export_1D(path=final_folder.joinpath("Si2O3.in"),
          lines=edit_line(path=initial_folder.joinpath("Si2O3.in"),
                          line=15,
                          replacement=f"read_restart\t\t\t{Path(initial_folder.name).joinpath('Si2O3_restart.restart')}"))
