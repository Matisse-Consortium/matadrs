from .tools import get_colormap

OPTIONS = {}

# NOTE: Averaging
OPTIONS["average.method"] = "avg_oifits"
OPTIONS["average.func"] = "robustmean"

# NOTE: Plot
OPTIONS["plot.colors.colormap"] = "seaborn-v0_8-colorblind"
OPTIONS["plot.colors.number"] = 10
OPTIONS["plot.color"] = get_colormap(OPTIONS["plot.colors.colormap"],
                                     OPTIONS["plot.colors.number"])
OPTIONS["plot.legend.fontsize"] = "small"
OPTIONS["plot.legend.location"] = "upper right"
OPTIONS["plot.linestyles"] = ["-", "--", "-.", ":"]
OPTIONS["plot.size"] = 700

# NOTE: Reduction
OPTIONS["reduce.cache.size"] = 1e4
