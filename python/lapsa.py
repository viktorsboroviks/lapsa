import pandas as pd
import vplot

DATA_COL_I = "i"
DATA_COL_T = "temperature"
DATA_COL_E = "energy"
COMPRESS_X_TO_N = 1000
SUBTITLE = "energy"
SUBTITLE_X = -0.03
SUBTITLE_Y = -0.07


def energy_temperature_subplot(
    data_csv,
    col,
    row,
    highlight_x=None,
    data_col_i=DATA_COL_I,
    data_col_t=DATA_COL_T,
    data_col_e=DATA_COL_E,
    compress_x_to_n=COMPRESS_X_TO_N,
    subtitle=SUBTITLE,
    subtitle_x=SUBTITLE_X,
    subtitle_y=SUBTITLE_Y,
):
    # read and compress energy table to save processing time
    try:
        table = pd.read_csv(data_csv)
        if table.index.size > compress_x_to_n:
            x_step = int(table.index.size / compress_x_to_n)
        else:
            x_step = 1
        compressed_table = table.iloc[::x_step]
    except FileNotFoundError:
        raise FileNotFoundError(f"plot_account_data skipped: {data_csv} not found")

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
        subtitle_text=subtitle,
        subtitle_x=subtitle_x,
        subtitle_y=subtitle_y,
        traces=[
            vplot.Step(
                x=compressed_table[data_col_i],
                y=compressed_table[data_col_t],
                color=vplot.Color.RED,
                name="temperature",
                showlegend=True,
            ),
            vplot.Step(
                x=compressed_table[data_col_i],
                y=compressed_table[data_col_e],
                color=vplot.Color.BLUE,
                name="energy",
                showlegend=True,
            ),
        ],
        lines=lines,
    )
