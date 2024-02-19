from pathlib import Path

from matadrs.utils.tools import query_archive


query_archive("MbS", instrument="matisse",
              target="hd142527", nights=[""],
              save_dir=Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142527/raw"))
