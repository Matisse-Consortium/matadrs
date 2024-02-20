from pathlib import Path

from matadrs.utils.query import query_archive


query_archive("MbS", instrument="matisse",
              target="hd142527", nights=["2021-03-08"], store_password=True,
              save_dir=Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142527/raw"))
