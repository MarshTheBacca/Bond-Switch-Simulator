from pathlib import Path

from utils import NetMCOutputData

cwd = Path(__file__).parent

data = NetMCOutputData(cwd.parent.joinpath("run", "output_files", "bss_stats.csv"))
data.summarise()
data.plot_energy()
