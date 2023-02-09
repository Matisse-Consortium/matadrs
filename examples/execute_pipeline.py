from pathlib import Path

from matadrs import reduction_pipeline


data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd163296")
observation_dirs = ["2019-03-22"]
raw_dirs = list(map(lambda x: data_dir / "raw" / x, observation_dirs))
product_dirs = list(map(lambda x: Path(str(x).replace("raw", "test")), raw_dirs))
reduction_pipeline(raw_dirs, product_dirs, overwrite=True, do_reduce=True,
                   do_calibrate=False, do_average=False, do_merge=False)
