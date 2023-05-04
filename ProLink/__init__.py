#!/usr/bin/env python3

import configparser
import os


# get path to ProLink folder
ProLink_path = os.path.dirname(os.path.realpath(__file__))

# read default parameters from configuration file (as a plain dictionary)
config = configparser.ConfigParser()
config.read(os.path.join(ProLink_path, 'parameters.cfg'))
parameters_default = {}
for section in config.sections():
    parameters_default.update(dict(config[section]))

# main function and parameters for better import and usage
from .prolink import pro_link
