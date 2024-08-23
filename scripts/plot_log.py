"""
Plot log.
"""

import argparse
import json

import vplot
import lapsa


CONFIG_SECTION = "plot_log"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--config", help="path to a .json config file")
args = parser.parse_args()

with open(args.config, encoding="UTF-8") as f:
    config_json = json.load(f)

calc_e_mean_n = config_json[CONFIG_SECTION]["calc_e_mean_n"]
calc_e_sd_n = config_json[CONFIG_SECTION]["calc_e_sd_n"]
log_pd = lapsa.log_to_pd(
    config_json[CONFIG_SECTION]["csv_path"],
    calc_e_mean_n=config_json[CONFIG_SECTION]["calc_e_mean_n"],
    calc_e_sd_n=config_json[CONFIG_SECTION]["calc_e_sd_n"],
)

subplots = []
subplots += [
    lapsa.subplot_e(
        log_pd,
        col=1,
        row=2,
        log_y=True,
        subtitle_text="e",
        calc_e_mean_n=calc_e_mean_n,
        calc_e_sd_n=calc_e_sd_n,
    )
]
subplots += [lapsa.subplot_do_cool(log_pd, col=1, row=3, subtitle_text="t-")]
subplots += [lapsa.subplot_t(log_pd, col=1, row=1, log_y=True, subtitle_text="t")]
subplots += [lapsa.subplot_v(log_pd, col=1, row=4, subtitle_text="v")]

vplot.PlotlyPlot(
    font_size=8,
    height=config_json[CONFIG_SECTION]["height"],
    width=config_json[CONFIG_SECTION]["width"],
    row_ratios=[2, 2, 1, 1],
    title_text="max_vector",
    subplots=subplots,
).to_file(
    filename=config_json[CONFIG_SECTION]["output_path"],
    scale=config_json[CONFIG_SECTION]["scale"],
)
