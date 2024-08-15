import pandas as pd
import vplot

DATA_COL_I = "i"
DATA_COL_T = "temperature"
DATA_COL_E = "energy"
COMPRESS_X_TO_N = 1000
SUBTITLE_TEXT = "energy"
SUBTITLE_X = -0.03
SUBTITLE_Y = -0.07


# read and compress energy table to save processing time
def _read_table(
    data_csv,
    compress_x_to_n,
) -> pd.DataFrame:
    try:
        table = pd.read_csv(data_csv)
        if table.index.size > compress_x_to_n:
            x_step = int(table.index.size / compress_x_to_n)
        else:
            x_step = 1
        return table.iloc[::x_step]
    except FileNotFoundError:
        raise FileNotFoundError(f"plot_account_data skipped: {data_csv} not found")


def energy_temperature_subplot(
    data_csv,
    col,
    row,
    log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_t=DATA_COL_T,
    data_col_e=DATA_COL_E,
    compress_x_to_n=COMPRESS_X_TO_N,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    table = _read_table(data_csv, compress_x_to_n)

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
        log_y=log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Step(
                x=table[data_col_i],
                y=table[data_col_t],
                color=vplot.Color.RED,
                name="temperature",
                showlegend=True,
            ),
            vplot.Step(
                x=table[data_col_i],
                y=table[data_col_e],
                color=vplot.Color.BLUE,
                name="energy",
                showlegend=True,
            ),
        ],
        lines=lines,
    )


def energy_subplot(
    data_csv,
    col,
    row,
    log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_e=DATA_COL_E,
    compress_x_to_n=COMPRESS_X_TO_N,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    table = _read_table(data_csv, compress_x_to_n)

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
        log_y=log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Step(
                x=table[data_col_i],
                y=table[data_col_e],
                color=vplot.Color.BLUE,
                name="energy",
                showlegend=True,
            ),
        ],
        lines=lines,
    )


def temperature_subplot(
    data_csv,
    col,
    row,
    log_y=False,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_t=DATA_COL_T,
    compress_x_to_n=COMPRESS_X_TO_N,
    subtitle_text=SUBTITLE_TEXT,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    table = _read_table(data_csv, compress_x_to_n)

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
        log_y=log_y,
        subtitle_text=subtitle_text,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Step(
                x=table[data_col_i],
                y=table[data_col_t],
                color=vplot.Color.RED,
                name="temperature",
                showlegend=True,
            ),
        ],
        lines=lines,
    )
