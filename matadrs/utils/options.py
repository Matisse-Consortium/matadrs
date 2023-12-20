OPTIONS = {}

# NOTE: Averaging
OPTIONS["average.method"] = "avg_oifits"
OPTIONS["average.func"] = "robustmean"

# NOTE: Plot
# TODO: Maybe use matplotlibcolors directly?
OPTIONS["plot.colors"] = ["b", "g", "r", "c", "m", "y"]
OPTIONS["plot.linestyles"] = ["-", "--", "-.", ":"]

# NOTE: Reduction
OPTIONS["reduce.cache.size"] = 1e4
