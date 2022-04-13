# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 17:42:02 2022

@author: eorland
"""

import argparse
import setup_functions
import ee_connect

parser = argparse.ArgumentParser(description='Handles directories for user.')
parser.add_argument('--filepath', type=str, required=True,
                   help='Full file path to project folder.')
parser.add_argument('-ap', '--archive_path', help='Location of output files\
                                                   including FIRMS data and asset downloads.\
                                                   Optional with default to the main directory if unspecified.')
parser.add_argument('-mp', '--model_output_path', help='Location of geojson output files from model.\
                                                        Optional with default to the main directory if unspecified')
parser.add_argument('--gee_username', type=str, required=True,
                   help='Case sensitive username for GEE account.')
# add specific path for model geojson output
args = parser.parse_args()
main_path = args.filepath

if args.archive_path:
    archive_path = args.archive_path
else: 
    archive_path = main_path

if args.model_output_path:
    model_output_path = args.model_output_path
else: 
    model_output_path = main_path

username = args.gee_username

ee_connect.connect()
setup_functions.directory_setup(archive_path, model_output_path)
setup_functions.make_records(main_path)
setup_functions.gee_folder_setup(username)