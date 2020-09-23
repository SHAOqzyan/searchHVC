from searchHVC  import   seachhvc

import os
import sys
from astropy.table import Table,vstack
from myPYTHON import *
from matplotlib.colors import LogNorm

doFITS=myFITS()
doHVC=seachhvc()


dataPath="/home/qzyan/WORK/projects/searchHVC/data/G210/"

if 0:
    doHVC.searchCloud(  "2100-035U.fits",path=dataPath )

if 1:

    allCOFITS=doHVC.getAllCO12FITS(dataPath)

    for eachCO in allCOFITS:
        print eachCO
        doHVC.searchCloud(eachCO,outPath=doHVC.tmpPath )

if 0:
    comBinTB=doHVC.getAllTBFiles(doHVC.tmpPath )

    #print len( comBinTB  )

if 0:
    doHVC.checkTB("combinedTB.fit")

if 0:

    TB=Table.read("combinedTB.fit")


    part1TB= doFITS.selectTBByColRange(TB,"v_cen", minV=90 )

    doHVC.drawCloudSpectra(part1TB ,dataPath)


    part2TB= doFITS.selectTBByColRange(TB,"v_cen", maxV=-200 )
    doHVC.drawCloudSpectra(part2TB ,dataPath)
