#!/usr/bin/env python
# coding: utf-8

# In[295]:

import urllib.request
import tempfile
import time
import os
import glob
import shutil
import requests, bs4, re, time
from selenium import webdriver
import numpy as np
import argparse
import sys
from selenium.webdriver.support.ui import WebDriverWait


##### FUNCTIONS #######
def str2bool(v):
    if isinstance(v, bool):
        return(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return(True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return(False)
    else:
        raise(argparse.ArgumentTypeError('Boolean value expected.'))

def pdbGenerator(SCOPid):
    
    sleepTime = 1 
    clickSleep = sleepTime*0.1
    exitCode = False 
    
    while exitCode == False: 
        driver = webdriver.Chrome()
        print('Waiting '+str(sleepTime)+' seconds for Web Page to load...')    
        driver.get('http://scop.mrc-lmb.cam.ac.uk/term/'+str(SCOPID))
        time.sleep(sleepTime)
        driver.implicitly_wait(sleepTime)
        like = driver.find_elements_by_class_name('icon')
        for x in range(0,len(like)):
            driver.implicitly_wait(sleepTime)
            if like[x].is_displayed():
                time.sleep(clickSleep)
                like[x].click()
        print('Waiting another '+str(sleepTime)+' seconds for Web Page to acquire source code...')
        time.sleep(sleepTime)
        jspagedata = bs4.BeautifulSoup(driver.page_source, 'html.parser')
        driver.close()

        repStructures = jspagedata.findAll('div',class_='grandchildrenlist')
        AAchainList = jspagedata.findAll('div', class_='td-protein-region')
        AAchainList = AAchainList[1::2]
        entries = len(repStructures)
        preprocessed_pdbName = []

        for i in range(entries):
            preprocessed_pdbName.extend(repStructures[i].findAll('div'))

        pdbName = [str(i)[5:9] for i in preprocessed_pdbName]
        pdbChain = [str(i)[10] for i in preprocessed_pdbName]
        
        #### Error Check for Downloads
        totalDomains = jspagedata.findAll('span', class_='title-note')
        totalDomains = np.int(str(totalDomains).split(' ')[-3])

        totalStructures = jspagedata.findAll('div', class_='descendants')
        totalStructures = np.sum([ np.int(str(numStruc).split(' ')[-2]) for numStruc in totalStructures])

        totDownloads = totalDomains+totalStructures
        
        ##################
    
        dat = jspagedata.find_all('a',href=True)
        preDomain = []
        for element in dat:
            if str(element).find('View in RCSB PDB') != -1:
                preDomain.append(element)

        representDomain = [str(pdb)[63:67] for pdb in preDomain]
        representDomainChain = [str(AA)[33:-6]+'_'+str(AA)[31:32] for AA in AAchainList]


        outputPDB = pdbName
        outputPDB.extend(representDomain)
        pdbChain.extend(representDomainChain)
        
        consistencyVal = len(outputPDB)

        if totDownloads == consistencyVal:
            print('Total of '+str(totDownloads)+' structures successfully scraped from SCOP.')
            exitCode = True
        elif sleepTime > 10:
            print('Wait time is too long. Will only export '+str(consistencyVal)+' out of '+str(totDownloads)+' potential structures.')
            exitCode = True
        else: 
            print('Total scraped structures ('+str(consistencyVal)+') inconsistent with total available ('+str(totDownloads)+').')
            print('Increasing the web page loading time to '+str(sleepTime+1)+' seconds and trying again.')
            sleepTime += 1 
            clickSleep = sleepTime*0.1
            exitCode = False 
    return(outputPDB,pdbChain)

def pdbDownloader(pdbList,pdbChain,SCOPid):
    size = len(pdbList)
    saveDir = os.getcwd()

    # You should change 'test' to your preferred folder.
    MYDIR = 'PDB_files_'+str(SCOPid)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print('created folder : ', MYDIR)

    else:
        print(MYDIR, 'folder already exists.')
    
    for i in range(size):
        nameIdentifier = str(pdbList[i]) + '_' + str(pdbChain[i])
        fileName = saveDir+ '/' +MYDIR+ '/' + nameIdentifier+'.pdb'
        httpAddress = 'https://files.rcsb.org/download/' + str(pdbList[i]) + '.pdb'
        urllib.request.urlretrieve(httpAddress, fileName)
            
        b = 'Downloading File: ' +str(i+1) +'/' + str(size)+ '    Current Progress: ' + str(round(i*100/(size-1),2)) +' %'
        print (b, end="\r")                                                                      

def pdbSlicer(SCOPid):
    MYDIR = 'PDB_files_'+str(SCOPid)
    pdbFiles = glob.glob(MYDIR+'/*.pdb')
    count = 1
    maxCount = len(pdbFiles)
    RECORDS = ['DBREF','DBREF1','DBREF2','SEQADV','SEQRES','MODRES','HET','HELIX','SHEET','SSBOND','LINK','CISPEP','SITE','ATOM','ANISOU','TER','HETATM']
    CHAINID = [[12], #DBREF
               [12], #DBREF1
               [12], #DBREF2
               [16], #SEQADV
               [11], #SEQRES
               [16], #MODRES
               [12], #HET
               [19,31], #HELIX
               [21,32], #SHEET
               [15,29], #SSBOND
               [21,51], #LINK
               [15,29], #CISPEP
               [22,33,44,55], #SITE
               [21], #ATOM
               [21], #ANISOU
               [21], #TER
               [21] #HETATM
              ]
    SEQID = [
            [[14,17],[20,23],[55,59],[62,66]], #DBREF
            [[14,17],[20,23]], #DBREF1
            [[45,54],[57,66]], #DBREF2
            [[18,21],[43,47]], #SEQADV
            None, #SEQRES
            [[18,21]], #MODRES
            [[13,16]], #HET
            [[21,24],[33,36]], #HELIX
            [[22,25],[33,36]], #SHEET
            [[17,20],[31,34]], #SSBOND
            [[22,25],[52,55]], #LINK
            [[17,20],[31,34]], #CISPEP
            None, #SITE
            [[22,25]], #ATOM
            [[22,25]], #ANISOU
            [[22,25]], #TER
            [[22,25]] #HETATM
            ]

         
    for pdbFile in pdbFiles:
        print('Processing on '+str(count)+' out of '+str(maxCount)+ ' files: '+str(pdbFile.split('\\')[-1]))
    
        logicTest = pdbFile.split('\\')[1][5:-4] #Takes the suffix of the file to determine if it is a chain or residue segment needed.
        if logicTest.isalpha() == True: #To separate the chains
            chain = pdbFile.split('_')[-1][0] #Chain of interest

            with open(pdbFile, "r") as f:
                lines = f.readlines()
            with open(pdbFile, "w") as f:
                for line in lines:
                    lineTest=line.split()[0] #First element of the PDB line
                    

                    if (lineTest in RECORDS) == False: #if not a line that includes a chain index just write to file.
                        f.write(line)
                    else:
                        lineHeaderIndex = RECORDS.index(lineTest) #Index of the line Header
                        
                        chainLineTest = [str(line[i]) for i in CHAINID[lineHeaderIndex]]
                        if chain in chainLineTest:
                            f.write(line)
                            
        else: # To separate the residues -- WORK IN PROGRESS
            residue = (pdbFile.split('_')[-2]).split('-')
            startResidue = np.int(residue[0])
            endResidue = np.int(residue[1])
            residueRange = np.arange(startResidue,endResidue+1,dtype=int)
            
            chain = str((pdbFile.split('_')[-1]).split('-')[0])[0]
            

            with open(pdbFile, "r") as f:
                lines = f.readlines()
            with open(pdbFile, "w") as f:
                for line in lines:
                    lineTest=line.split()[0]
                    chainTest = False

                    if ((lineTest in RECORDS) == False):
                        f.write(line)
                    elif RECORDS.index(lineTest) in [4,12]:
                        f.write(line)
                    elif RECORDS.index(lineTest) in [0,1,2,3,5,6,7,8,9,10,11,13,14,15,16]:
                        lineHeaderIndex = RECORDS.index(lineTest) #Index of the line Header
                        
                        chainLineTest = [str(line[i]) for i in CHAINID[lineHeaderIndex]]
                        if chain in chainLineTest:
                            chainTest = True 

                        startSeqLineTest = [str(line[indices[0]:indices[1]+1]).strip() for indices in SEQID[lineHeaderIndex]]

                        if (lineHeaderIndex == 0) and chainTest: #DBREF
                            if startSeqLineTest[0].isnumeric() and startSeqLineTest[1].isnumeric():
                                recordResidueRange_first = np.arange(np.int(startSeqLineTest[0]),np.int(startSeqLineTest[1])+1,dtype=int)
                                if (np.intersect1d(residueRange,recordResidueRange_first).size != 0):
                                    f.write(line)
                            if startSeqLineTest[2].isnumeric() and startSeqLineTest[3].isnumeric():
                                recordResidueRange_second = np.arange(np.int(startSeqLineTest[2]),np.int(startSeqLineTest[3])+1,dtype=int)
                                if (np.intersect1d(residueRange,recordResidueRange_second).size != 0):
                                    f.write(line)
                        elif (lineHeaderIndex in [1,2,3,7,8,9,10,11]) and chainTest:
                            if startSeqLineTest[0].isnumeric() and startSeqLineTest[1].isnumeric():
                                recordResidueRange = np.arange(np.int(startSeqLineTest[0]),np.int(startSeqLineTest[1])+1,dtype=int)
                                if np.intersect1d(residueRange,recordResidueRange).size != 0:
                                    f.write(line)
                        elif chainTest:
                            if startSeqLineTest[0].isnumeric():
                                if np.int(startSeqLineTest[0]) in residueRange:
                                    f.write(line)
        count+=1
        
##########3
# In[326]:


parser = argparse.ArgumentParser(description='Downloads and zips all the PDBs associated with the SCOP ID.')
parser.add_argument('SCOPid', help='Enter the SCOP ID', type=int)
parser.add_argument("-splice", type=str2bool, nargs='?', const=True, default=False, help='Excise the chains/residues of interest from the PDB files. Work in Progress.')
parser.add_argument("-download", type=str2bool, nargs='?', const=True, default=True, help='Scrape and download the PDBs listed for the SCOPid.')
parser.add_argument("-archive", type=str2bool, nargs='?', const=True, default=True, help='Archive the folder into a zip file.')
args = parser.parse_args()

SCOPID = args.SCOPid

if str(SCOPID).isnumeric() != True:
    print('SCOP ID is not valid')
    sys.exit()

if args.download == True:
    test0, test1 = pdbGenerator(SCOPID)
    pdbDownloader(test0,test1,SCOPID)
    b = 'Downloading File: ' +str(len(test0)) +'/' + str(len(test0))+ '    Current Progress: ' + str(round(100,2)) +' %'
    print (b)
    print('Downloading Files Complete.')
    

if args.splice == True:
    pdbSlicer(SCOPID)
    print('Processing Complete.')
    
if args.archive == True:
    shutil.make_archive('PDB_files_'+str(SCOPID),'zip','PDB_files_'+str(SCOPID))
    print('Files are now archived in '+'PDB_files_'+str(SCOPID)+'.zip')

print('All processes are now complete.')

