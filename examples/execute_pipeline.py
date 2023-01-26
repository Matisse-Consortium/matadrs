from pathlib import Path

from matadrs import reduction_pipeline


data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666")
observation_dirs = ["UTs/2019-05-14", "UTs/2022-04-20", "UTs/2022-04-22"]
raw_dirs = list(map(lambda x: data_dir / "raw" / x, observation_dirs))
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "products")), raw_dirs))
reduction_pipeline(raw_dirs[1], product_dirs[1], overwrite=True, do_reduce=True,
                   do_calibrate=False, do_average=False, do_merge=False)
