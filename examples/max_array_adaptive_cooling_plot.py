"""
Plot max array - adaptive cooling.
"""

import argparse
import json
import vplot
import lapsa


DEFAULT_CONFIG_SECTION = "plot"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--config", help="path to a .json config file")
parser.add_argument(
    "--config-section",
    help="name of the section in config",
    default=DEFAULT_CONFIG_SECTION,
    required=False,
)
args = parser.parse_args()

with open(args.config, encoding="UTF-8") as f:
    config_json = json.load(f)
if args.config_section:
    config_json = config_json[args.config_section]

data_pd = lapsa.log_to_pd(config_json["data_csv"])
subplots = []
subplots += [
    lapsa.subplot_t(
        data_pd=data_pd,
        col=1,
        row=1,
        plot_log_y=True,
        data_col_i="run_i",
        subtitle_text=None,
    )
]

subplots += [
    lapsa.subplot_e(
        data_pd=data_pd,
        col=1,
        row=2,
        plot_log_y=True,
        data_col_i="run_i",
        subtitle_text=None,
    )
]

vplot.PlotlyPlot(
    title_text="max array (adaptive cooling)",
    subplots=subplots,
).to_file(config_json["output_file"])
