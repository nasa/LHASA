# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:51:14 2022

@author: eorland
"""

import argparse
import firms_basin_join

parser = argparse.ArgumentParser(description='Handles directories for user.')
parser.add_argument('--filepath', type=str, required=True,
                   help='Full file path of project folder')
parser.add_argument('--firms_path', type=str, required=True,
                   help='Full path to FIRMS data.')
parser.add_argument('--gee_username', type=str, required=True,
                   help='Case sensitive username for GEE account.')


args = parser.parse_args()
main_path = args.filepath
firms_path = args.firms_path
username = args.gee_username
print('join.py args:', args)
print('Calling firms_basin_join.join()...')
basins, run_dates, detects = firms_basin_join.join(main_path, firms_path)