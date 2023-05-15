from pathlib import Path

from matadrs.utils.plot import Plotter

file = Path("/Users/scheuck/Data/reduced_data/")

plotter = Plotter(file)
plotter.add_uv().add_vis().plot()
