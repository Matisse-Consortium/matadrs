from types import SimpleNamespace


# NOTE: Averaging
average = SimpleNamespace(method="avg_oifits",
                          func="robustmean")

# NOTE: Plot
color = SimpleNamespace(colormap="seaborn-v0_8-colorblind",
                        number=100)
legend = SimpleNamespace(fontsize="small",
                         location="upper right")

plot = SimpleNamespace(
        color=color, dpi=300,
        legend=legend,
        linestyles=["-", "--", "-.", ":"],
        size=700)

OPTIONS = SimpleNamespace(
        average=average,
        color=color,
        legend=legend,
        plot=plot)

