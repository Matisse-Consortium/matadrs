from pathlib import Path

from matadrs import reduction_pipeline


# Specify the path to the directory containing the data
data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666/")
# observation_dirs = [data_dir / "raw/UTs/2022-04-21",
                    # data_dir / "raw/UTs/2022-04-23"]

observation_dirs = [data_dir / "raw/UTs/2019-05-13"]

# Specify the raw-directory, containing the raw data
# raw_dirs = list(map(lambda x: data_dir / "raw" / x, observation_dirs))

# Specify the product-directory, to contain the product data/that contains reduced,
# calibrated or averaged data, to be further processed
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "product")),
                        observation_dirs))

# Call the reduction_pipeline
reduction_pipeline(observation_dirs, product_dirs, overwrite=True, do_reduce=False,
                   do_calibrate=True, do_average=False, do_merge=False)
