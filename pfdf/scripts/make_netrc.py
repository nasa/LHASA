"""
The following script helps set up a .netrc file to access HLS data from LPDAAC 
via AWS. This code is based on a deleted page on the NASA EarthData Cloud 
Cookbook:

https://nasa-openscapes.github.io/earthdata-cloud-cookbook/
"""

from netrc import netrc
from platform import system
from getpass import getpass
import os

urs = "urs.earthdata.nasa.gov"  # Earthdata URL endpoint for authentication
prompts = [
    "Enter NASA Earthdata Login Username: ",
    "Enter NASA Earthdata Login Password: ",
]

# Determine the OS (Windows machines usually use an '_netrc' file)
netrc_name = "_netrc" if system() == "Windows" else ".netrc"
file_path = os.path.join(os.path.expanduser("~"), netrc_name)

# Determine if netrc file exists, and if so, if it includes NASA Earthdata Login Credentials
try:
    netrc(file_path).authenticators(urs)[0]
    print(f"{file_path} already contains valid NASA Earthdata Login.")
except:
    text = (
        f"machine {urs}\n"
        f"  login {getpass(prompt=prompts[0])}\n"
        f"  password {getpass(prompt=prompts[1])}\n"
    )
    print(f"Adding NASA Earthdata Login Credentials to {file_path}")
    with open(file_path, mode="a") as f:
        os.chmod(file_path, 0o600)
        f.write(text)
