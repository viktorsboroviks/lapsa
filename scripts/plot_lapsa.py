import argparse
import json
import pandas as pd
import vplot

# pylint: disable=duplicate-code
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--config", help="path to a .json config file")
args = parser.parse_args()

with open(args.config, encoding="UTF-8") as f:
    config_json = json.load(f)
# pylint: enable=duplicate-code

output_file = config_json["plot_lapsa"]["output_file_name"]
log_df = pd.read_csv(config_json["lapsa"]["log_file_name"])

subplots = []
subplots += [
    vplot.Subplot(
        col=1,
        row=1,
        log_y=True,
        x_title="run_i",
        y_title="t",
        traces=[
            vplot.Scatter(
                x=log_df["run_i"],
                y=log_df["temperature"],
                color=vplot.Color.RED,
                mode="lines",
                name="t",
            )
        ],
    )
]

subplots += [
    vplot.Subplot(
        col=1,
        row=2,
        log_y=True,
        x_title="run_i",
        y_title="e",
        traces=[
            vplot.Scatter(
                x=log_df["run_i"],
                y=log_df["energy"],
                color=vplot.Color.BLUE,
                mode="lines",
                name="e",
            ),
        ],
    )
]

subplots += [
    vplot.Subplot(
        col=1,
        row=3,
        log_y=True,
        x_title="run_i",
        y_title="v",
        traces=[
            vplot.Scatter(
                x=log_df["run_i"],
                y=log_df["value"],
                color=vplot.Color.BLUE,
                mode="lines",
                name="v",
            ),
        ],
    )
]

vplot.PlotlyPlot(
    subplots=subplots,
).to_file(output_file)
