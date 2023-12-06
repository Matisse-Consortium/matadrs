from pathlib import Path

from matadrs.utils.plot import Plotter, plot_data_quality

# file = Path("/Users/scheuck/Data/calibrated_data/lband/test/TARGET_FLUX_CALT_0004.fits")

# plotter = Plotter(file)
# plotter.add_uv().add_vis().plot(error=True)


path = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/"
            "data/matisse/GTO/hd142666/product/UTs/")
for directory in list(path.glob("20*")):
    reduced_dir = directory / "reduced" / "incoherent" / "lband"
    plot_data_quality(reduced_dir, reduced_dir / "data_quality")
