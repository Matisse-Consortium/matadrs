from pathlib import Path

from matadrs.utils.plot import Plotter

file = Path("/Users/scheuck/Data/calibrated_data/lband/test/TARGET_FLUX_CALT_0004.fits")

plotter = Plotter(file)
plotter.add_uv().add_vis().plot(error=True)
