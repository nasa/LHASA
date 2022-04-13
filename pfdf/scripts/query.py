# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 10:45:31 2022

@author: eorland
"""

import argparse
import ee_connect
from query_and_download import query_and_download

parser = argparse.ArgumentParser(description='Handles directories for user.')
parser.add_argument('--filepath', type=str, required=True,
                   help='Full file path of project folder')
parser.add_argument('-ap', '--asset_path',help='Path to folder containing\
                                                the directory "asset_downloads".')
parser.add_argument('--gee_username', type=str, required=True,
                   help='Case sensitive username for GEE account.')

args = parser.parse_args()
main_path = args.filepath

if args.asset_path:
    asset_path = args.asset_path
else: 
    asset_path = main_path

username = args.gee_username

ee_connect.connect()
query_and_download(main_path,asset_path,username)