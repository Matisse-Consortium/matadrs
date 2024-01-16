OPTIONS = {}

# NOTE: Averaging
OPTIONS["average.method"] = "avg_oifits"
OPTIONS["average.func"] = "robustmean"

# NOTE: Plot
OPTIONS["plot.dpi"] = 300
OPTIONS["plot.color.colormap"] = "seaborn-v0_8-colorblind"
OPTIONS["plot.color.number"] = 100
OPTIONS["plot.legend.fontsize"] = "small"
OPTIONS["plot.legend.location"] = "upper right"
OPTIONS["plot.linestyles"] = ["-", "--", "-.", ":"]
OPTIONS["plot.size"] = 700

# NOTE: Reduction
OPTIONS["reduce.cache.size"] = 1e4
