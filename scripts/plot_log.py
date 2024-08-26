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
row_ratios = []
if config_json[CONFIG_SECTION]["plot_t"]:
    subplots += [
        lapsa.subplot_t(log_pd, col=1, row=1, plot_log_y=True, subtitle_text="t")
    ]
    row_ratios += [2]
if config_json[CONFIG_SECTION]["plot_do_cool"]:
    subplots += [lapsa.subplot_do_cool(log_pd, col=1, row=2, subtitle_text="t-")]
    row_ratios += [1]
if config_json[CONFIG_SECTION]["plot_e"]:
    subplots += [
        lapsa.subplot_e(
            log_pd,
            col=1,
            row=3,
            plot_log_y=True,
            plot_mean=bool(calc_e_mean_n),
            plot_sd=bool(calc_e_sd_n),
            subtitle_text="E",
        )
    ]
    row_ratios += [2]
if config_json[CONFIG_SECTION]["plot_de"]:
    subplots += [
        lapsa.subplot_de(
            log_pd,
            col=1,
            row=4,
            plot_log_y=False,
            plot_mean=bool(calc_e_mean_n),
            plot_sd=bool(calc_e_sd_n),
            subtitle_text="dE",
        )
    ]
    row_ratios += [2]
if config_json[CONFIG_SECTION]["plot_v"]:
    subplots += [lapsa.subplot_v(log_pd, col=1, row=5, subtitle_text="v")]
    row_ratios += [1]

vplot.PlotlyPlot(
    font_size=8,
    height=config_json[CONFIG_SECTION]["height"],
    width=config_json[CONFIG_SECTION]["width"],
    title_text="max_vector",
    subplots=subplots,
    row_ratios=row_ratios,
).to_file(
    filename=config_json[CONFIG_SECTION]["output_path"],
    scale=config_json[CONFIG_SECTION]["scale"],
)
