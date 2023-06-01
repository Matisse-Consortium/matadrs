from pathlib import Path

from matadrs.utils.readout import ReadoutFits

file = Path("/Users/scheuck/Data/reduced_data/hd142666/midi/hd142666_2004-06-08T04:19:12_MIDI.fits")

readout = ReadoutFits(file)
print(readout.oi_wl)
print(readout.oi_vis)
print(readout.oi_vis.columns)
