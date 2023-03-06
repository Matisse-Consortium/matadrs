from pathlib import Path

from matadrs.utils.plot import Plotter

file = Path("/Users/scheuck/Data/reduced_data/hd142666/midi/hd142666_2004-06-08T04:19:12_MIDI.fits")

plotter = Plotter(file)
plotter.add_uv().add_vis().plot()
