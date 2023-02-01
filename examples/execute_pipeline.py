from pathlib import Path

from matadrs import reduction_pipeline


data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666")
observation_dirs = ["UTs/2022-04-23"]
raw_dirs = list(map(lambda x: data_dir / "raw" / x, observation_dirs))
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "test")), raw_dirs))
reduction_pipeline(raw_dirs, product_dirs, overwrite=True, do_reduce=True,
                   do_calibrate=True, do_average=True, do_merge=True)
