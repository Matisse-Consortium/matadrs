from pathlib import Path

from matadrs import matadrs_pipeline
from matadrs.utils.options import OPTIONS


# Specify the path to the directory containing the data
data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142527")

# Specify the raw-directory, containing the raw data
# raw_dirs = [data_dir / "raw" / "UTS"]
observation_dirs = data_dir / "raw" / "test"
raw_dirs = sorted(observation_dirs.glob("*"), key=lambda x: x.name[-8:])

# Specify the product-directory, to contain the product data/that contains reduced,
# calibrated or averaged data, to be further processed
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "product")), raw_dirs))

# Specify averaging pipeline
OPTIONS["average.method"] = "mat_tools"

# Call the reduction_pipeline
matadrs_pipeline(raw_dirs, product_dirs, overwrite=True, band="nband",
                 do_reduce=True, do_calibrate=False,
                 do_average=False, do_merge=False, ncores=6)
