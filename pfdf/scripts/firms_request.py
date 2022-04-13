# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:27:57 2022

@author: eorland
"""

import os
import requests
import bs4

def download_firms(main_path,firms_path,show_progress=False,overwrite=False):
    '''
    Downloads recent FIRMS data from the following link: 
    https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/
    
    Takes the project directory as an input and saves all files in the 
    corresponding directory firms/all_firms'''
    
    token_path = os.path.join(main_path, 'ref_data','token.txt')
    
    with open(token_path) as f:
        token = f.read().strip()
    
   # token = 'this is my stupid token twueouws aSabsagjskjaS'
    print('using token:', token)
    
    # build path to download directory
    os.chdir(firms_path)
    
    base = 'https://nrt3.modaps.eosdis.nasa.gov' # host url
    # urls for each VIIRS sensor
    soumi_url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/suomi-npp-viirs-c2/Global'
    noaa_url = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/noaa-20-viirs-c2/Global'
    url_list = [soumi_url,noaa_url]
    
    for url in url_list:
        r = requests.get(url)#,headers =
                        #{'Authorization': 'Bearer {}'.format(token)})
        
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
            r = requests.get(link,timeout=5,headers =
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
            r = requests.get(link,timeout=5,headers =
                             {'Authorization': 'Bearer {}'.format(token)})
            if r.ok:
                with open(file_name, "w") as file:
                    if show_progress:
                        print("Downloading", file_name + '...')
                    file.write(r.text)
            else: 
                print('\nCould not download', file_name +
                      ' \nHTTP status code:', r.status_code)
                
        