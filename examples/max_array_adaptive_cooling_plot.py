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

with open(args.config) as f:
    config_json = json.load(f)
if args.config_section:
    config_json = config_json[args.config_section]

subplot = lapsa.energy_temperature_subplot(
    data_csv=config_json["data_csv"], col=1, row=1, data_col_i="state_i", subtitle=None
)

vplot.PlotlyPlot(title_text="max array (adaptive cooling)", subplots=[subplot]).to_file(
    config_json["output_file"]
)
