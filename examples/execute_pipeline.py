from pathlib import Path

from matadrs import matadrs_pipeline
from matadrs.utils.options import OPTIONS


# Specify the path to the directory containing the data
data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666")

# Specify the raw-directory, containing the raw data
observation_dirs = data_dir / "raw" / "UTs"
raw_dirs = sorted(observation_dirs.glob("20*"), key=lambda x: x.name[-8:])

# Specify the product-directory, to contain the product data/that contains reduced,
# calibrated or averaged data, to be further processed
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "product")), raw_dirs))

# Specify averaging pipeline
OPTIONS["average.method"] = "mat_tools"

# Call the reduction_pipeline
matadrs_pipeline(raw_dirs, product_dirs, overwrite=True,
                 do_reduce=False, do_calibrate=False,
                 do_average=True, do_merge=True, ncores=6)
