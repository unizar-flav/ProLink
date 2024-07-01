#!/usr/bin/env python3


import logging
import os

import yaml


__version__ = '0.2.0'

# set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[
        logging.StreamHandler()
        ]
    )

# get path to ProLink folder
ProLink_path = os.path.dirname(os.path.realpath(__file__))

# read default parameters from configuration file (as a plain dictionary)
with open(os.path.join(ProLink_path, 'parameters.yaml'), 'r') as f:
    parameters_yaml = yaml.safe_load(f)
    parameters_default = {key: value for section in parameters_yaml.values() for key, value in section.items()}

# main function and parameters for better import and usage
from .prolink import pro_link, pro_link_multiple
