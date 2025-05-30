"""
Lapsa plots.
"""

# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import dataclasses
import pandas as pd
import vplot


DATA_COL_I = "run_i"
DATA_COL_T = "temperature"
DATA_COL_E = "energy"
DATA_COL_V = "value"
DATA_COL_E_MEAN = "energy mean"
DATA_COL_E_SD = "energy sd"
DATA_COL_DE = "d energy"
DATA_COL_DE_MEAN = "d energy mean"
DATA_COL_DE_SD = "d energy sd"
DATA_COL_DE_SD_MEAN = "d energy sd mean"
DATA_COL_DO_COOL = "do cool"
SUBTITLE_TEXT = "energy"
SUBTITLE_X = -0.03
SUBTITLE_Y = -0.07


# read and compress energy table to save processing time
def log_to_pd(
    log_csv,
    calc_e_mean_n=None,
    calc_e_sd_n=None,
    compress_x_to_n=None,
    data_col_e=DATA_COL_E,
    data_col_e_mean=DATA_COL_E_MEAN,
    data_col_e_sd=DATA_COL_E_SD,
    data_col_de=DATA_COL_DE,
    data_col_de_mean=DATA_COL_DE_MEAN,
    data_col_de_sd=DATA_COL_DE_SD,
    data_col_t=DATA_COL_T,
    data_col_do_cool=DATA_COL_DO_COOL,
) -> pd.DataFrame:
    try:
        table = pd.read_csv(log_csv)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"plot_account_data aborted: {log_csv} not found"
        ) from e

    # calculate when temperature was decreased
    table[data_col_do_cool] = table[data_col_t] < table[data_col_t].shift(1)

    # de
    table[data_col_de] = table[data_col_e] - table[data_col_e].shift(1)

    # derive extra data
    # - mean
    # - standard deviation
    if calc_e_mean_n:
        table[data_col_e_mean] = table[data_col_e].rolling(window=calc_e_mean_n).mean()
        table[data_col_de_mean] = (
            table[data_col_de].rolling(window=calc_e_mean_n).mean()
        )
    if calc_e_sd_n:
        table[data_col_e_sd] = table[data_col_e].rolling(window=calc_e_sd_n).std()
        table[data_col_de_sd] = table[data_col_de].rolling(window=calc_e_sd_n).std()

    # compress the table to a smaller format to speed up rendering
    if compress_x_to_n and table.index.size > compress_x_to_n:
        x_step = int(table.index.size / compress_x_to_n)
    else:
        x_step = 1
    return table.iloc[::x_step]


def subplot_t(
    data_pd,
    col,
    row,
    plot_log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_t=DATA_COL_T,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    lines = None
    if highlight_x:
        lines = [
            vplot.Lines(
                x=[highlight_x],
                color=vplot.Color.RED,
                dash=vplot.Dash.SOLID,
            )
        ]

    return vplot.Subplot(
        col=col,
        row=row,
        log_y=plot_log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Step(
                x=data_pd[data_col_i],
                y=data_pd[data_col_t],
                color=vplot.Color.RED,
                name="temperature",
                showlegend=True,
            ),
        ],
        lines=lines,
    )


def subplot_e(
    data_pd,
    col,
    row,
    name="E",
    plot_data=True,
    plot_mean=False,
    plot_sd=False,
    plot_log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_e=DATA_COL_E,
    data_col_e_mean=DATA_COL_E_MEAN,
    data_col_e_sd=DATA_COL_E_SD,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    lines = None
    if highlight_x:
        lines = [
            vplot.Lines(
                x=[highlight_x],
                color=vplot.Color.RED,
                dash=vplot.Dash.SOLID,
            )
        ]

    traces = []

    if plot_data:
        traces += [
            vplot.Scatter(
                x=data_pd[data_col_i],
                y=data_pd[data_col_e],
                color=vplot.Color.BLUE,
                name=name,
                showlegend=True,
            )
        ]

    if plot_mean:
        traces += [
            vplot.Scatter(
                x=data_pd[data_col_i],
                y=data_pd[data_col_e_mean],
                color=vplot.Color.RED,
                dash=vplot.Dash.SOLID,
                name="mean",
                showlegend=True,
            )
        ]

    if plot_sd:

        @dataclasses.dataclass
        class SD:
            coef: float
            dash: vplot.Dash
            name: str

        sds = [
            SD(coef=+2, dash=vplot.Dash.DOT, name="+2sd"),
            SD(coef=+1, dash=vplot.Dash.SOLID, name="+1sd"),
            SD(coef=-1, dash=vplot.Dash.SOLID, name="-1sd"),
            SD(coef=-2, dash=vplot.Dash.DOT, name="-2sd"),
        ]

        for sd in sds:
            traces += [
                vplot.Scatter(
                    x=data_pd[data_col_i],
                    y=data_pd[data_col_e_mean] + sd.coef * data_pd[data_col_e_sd],
                    color=vplot.Color.ORANGE,
                    dash=sd.dash,
                    name=sd.name,
                    showlegend=True,
                )
            ]

    return vplot.Subplot(
        col=col,
        row=row,
        log_y=plot_log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=traces,
        lines=lines,
    )


def subplot_de(
    data_pd,
    col,
    row,
    name="dE",
    plot_data=True,
    plot_mean=False,
    plot_sd=False,
    plot_log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_e=DATA_COL_DE,
    data_col_e_mean=DATA_COL_DE_MEAN,
    data_col_e_sd=DATA_COL_DE_SD,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    return subplot_e(
        data_pd=data_pd,
        col=col,
        row=row,
        name=name,
        plot_data=plot_data,
        plot_mean=plot_mean,
        plot_sd=plot_sd,
        plot_log_y=plot_log_y,
        highlight_x=highlight_x,
        data_col_i=data_col_i,
        data_col_e=data_col_e,
        data_col_e_mean=data_col_e_mean,
        data_col_e_sd=data_col_e_sd,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
    )


def subplot_do_cool(
    data_pd,
    col,
    row,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_do_cool=DATA_COL_DO_COOL,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    lines = None
    if highlight_x:
        lines = [
            vplot.Lines(
                x=[highlight_x],
                color=vplot.Color.RED,
                dash=vplot.Dash.SOLID,
            )
        ]

    traces = []
    traces += [
        vplot.Step(
            x=data_pd[data_col_i],
            y=data_pd[data_col_do_cool],
            color=vplot.Color.BLUE,
            name="do cool",
            showlegend=True,
        )
    ]

    return vplot.Subplot(
        col=col,
        row=row,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=traces,
        lines=lines,
    )


def subplot_v(
    data_pd,
    col,
    row,
    plot_log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_v=DATA_COL_V,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    lines = None
    if highlight_x:
        lines = [
            vplot.Lines(
                x=[highlight_x],
                color=vplot.Color.RED,
                dash=vplot.Dash.SOLID,
            )
        ]

    return vplot.Subplot(
        col=col,
        row=row,
        log_y=plot_log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Scatter(
                x=data_pd[data_col_i],
                y=data_pd[data_col_v],
                color=vplot.Color.BLACK,
                name="value",
                showlegend=True,
            ),
        ],
        lines=lines,
    )
