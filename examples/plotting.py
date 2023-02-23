from matadrs.utils.plot import Plotter

data = "/Users/scheuck/Data/reduced_data/jozsef_reductions/hd142666/HD_142666_2022-04-21T07_18_22_N_TARGET_FINALCAL_INT.fits"

plotter = Plotter(data)
plotter.add_uv().add_vis().add_vis(corr_flux=True)
plotter.add_flux().add_cphases(unwrap=True).add_diff_phases(unwrap=True).plot()
