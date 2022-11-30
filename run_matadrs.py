from matadrs import matadrs_pipeline

DATA_DIR = "~/Data/hd163296/"
STEM_DIR = ""
TARGET_LST = ["20190323"]
TARGET_DIRS = [os.path.join("ATs", folder) for folder in TARGET_LST]
matadrs_pipeline(DATA_DIR, STEM_DIR, TARGET_DIRS,
                 do_reduce=False, do_calibrate=True,
                 do_average=False, do_merge=False)
