from utils import NetMCOutputData
from pathlib import Path

cwd = Path(__file__).parent

data = NetMCOutputData(cwd.parent.joinpath("testing","output_files","test_all_stats.csv"))
data.plot_ring_sizes()