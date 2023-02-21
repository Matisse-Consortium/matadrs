import matplotlib.pyplot as plt

from matadrs.utils.plot import Plotter

data = "/Users/scheuck/Data/reduced_data/jozsef_reductions/hd142666/HD_142666_2022-04-21T07_18_22_L_TARGET_FINALCAL_INT.fits"

plotter = Plotter(data)
fig, ax = plt.subplots(1, 1)
plotter.plot_uv(ax)
plt.show()
