# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:51:49 2022

@author: eorland
"""

import argparse
import firms_request

parser = argparse.ArgumentParser(description='Handles directories for user.')
parser.add_argument('--filepath', type=str, required=True,
                   help='Full file path of project folder')
parser.add_argument('-fp', '--firms_path', required=True, help='location of FIRMS files')

args = parser.parse_args()
main_path = args.filepath
firms_path = args.firms_path

firms_request.download_firms(main_path,firms_path,show_progress=True,overwrite=False)