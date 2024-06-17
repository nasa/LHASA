# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:27:57 2022

@author: eorland
"""

import os
import requests
import bs4
import argparse

def download_firms(home_path,firms_path,show_progress=False,overwrite=False):
    '''
    Downloads recent FIRMS data from the following link: 
    https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/
    
    Takes the project directory as an input and saves all files in the 
    corresponding directory firms/all_firms
    '''
    
    token_path = os.path.join(home_path, 'ref_data','token.txt')
    
    with open(token_path) as f:
        token = f.read().strip()
    
    os.chdir(firms_path)
    
    base_url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/'
    soumi_url = f'{base_url}suomi-npp-viirs-c2/Global'
    noaa_url = f'{base_url}noaa-20-viirs-c2/Global'
    url_list = [soumi_url, noaa_url]
    
    for url in url_list:
        r = requests.get(url)
        data = bs4.BeautifulSoup(r.text, "html.parser")
        
        # loop through all <a> tags - these are the file paths 
        for l in data.find_all("a")[1:-2]: # index [0] is NOT a file path
            link = l["href"] # base link + file link = full url 
            file_name = os.path.basename(link)
          
            if show_progress and os.path.exists(file_name):
                print(file_name, "already exists.")
                if overwrite:
                    print('Overwriting',file_name)
                    r = requests.get(link,headers =
                             {'Authorization': 'Bearer {}'.format(token)})
                    if r.ok:
                        with open(file_name, "w") as file:
                            file.write(r.text)
                    else: 
                        print('\nCould not download', file_name +
                              ' \nHTTP status code:', r.status_code)
                else: 
                    continue
          
            r = requests.get(link,headers =
                             {'Authorization': 'Bearer {}'.format(token)})
            if r.ok:
                with open(file_name, "w") as file:
                    if show_progress:
                        print("Downloading", file_name + '...')
                    file.write(r.text)
            else: 
                print('\nCould not download', file_name +
                      ' \nHTTP status code:', r.status_code)
        
        # last step, redownload the previous day's file
        # this is to ensure we get the final data for it vs provisional data
        for l in data.find_all("a")[-2:]: # indices are last two -- always downloaded
            link = l["href"] # base link + file link = full url 
            file_name = os.path.basename(link)
            r = requests.get(link,headers =
                             {'Authorization': 'Bearer {}'.format(token)})
            if r.ok:
                with open(file_name, "w") as file:
                    if show_progress:
                        print("Downloading", file_name + '...')
                    file.write(r.text)
            else: 
                print('\nCould not download', file_name +
                      ' \nHTTP status code:', r.status_code)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Handles directories for user.')
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="LHASA version 2.0.1",
    )
    parser.add_argument(
        "-p",
        "--path",
        default=os.getcwd(),
        help="file path to project folder.",
    )
    parser.add_argument(
        "-ap",
        "--archive_path",
        help=(
            "Location of output files including FIRMS data and asset downloads"
            ". Optional with default to the main directory if unspecified."
        ),
    )

    args = parser.parse_args()
    home_path = args.path
    if args.archive_path:
        firms_path = os.path.join(args.archive_path, "firms", "all_firms")
    else:
        firms_path = os.path.join(home_path, "firms", "all_firms")

    download_firms(home_path,firms_path,show_progress=True,overwrite=False)
