Needs a webdriver for your browser to work. 

If using chrome, go to 
https://chromedriver.chromium.org/downloads
to download a webdriver for your browser version.

The webdriver needs to match your browser version as well as the operating system being used. 

PLEASE PLACE THE WEBDRIVER IN THE SAME LOCATION AS THE CODE.

If you're running on a Unix-based operating system, you may need to update the permissions of chromedriver after downloading it in order to make it executable:

chmod +x chromedriver

HOW TO USE THIS

> python scop.py SCOPid [-h] [-splice [bool]] [-download [bool]] [-archive [bool]]

Downloads and zips all the PDBs associated with the SCOP ID.

positional arguments:
  SCOPid                Enter the SCOP ID

optional arguments:
  -h, --help            show this help message and exit
  -splice [SPLICE]      Excise the chains/residues of interest from the PDB files. This is still a work in progress, please double check                         your files and use at your own risk.
  -download [DOWNLOAD]  Scrape and download the PDBs listed for the SCOPid.
  -archive [ARCHIVE]    Archive the folder into a zip file.
