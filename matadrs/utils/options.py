from types import SimpleNamespace


# NOTE: Averaging
average = SimpleNamespace(
    method="mat_tools",
    func="robustmean"
    )

# NOTE: Plot
color = SimpleNamespace(
        colormap="seaborn-v0_8-colorblind",
        number=35
        )

legend = SimpleNamespace(
        fontsize="small",
        location="upper right"
        )

plot = SimpleNamespace(
        color=color, dpi=300,
        legend=legend,
        linestyles=["-", "--", "-.", ":"],
        size=700
        )

# NOTE: All settings
OPTIONS = SimpleNamespace(
        average=average,
        color=color,
        legend=legend,
        plot=plot
        )
