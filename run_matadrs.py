import os

from matadrs import matadrs_pipeline

DATA_DIR = "/Users/scheuck/Data/data_reduction/"
STEM_DIR = "hd163296/"
TARGET_DIRS = ["ATs/20190323"]
matadrs_pipeline(DATA_DIR, STEM_DIR, TARGET_DIRS,
                 do_reduce=True, do_calibrate=True,
                 do_average=False, do_merge=False)


