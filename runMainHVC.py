from searchHVC  import   seachhvc

import os
import sys
from astropy.table import Table,vstack
from myPYTHON import *
from matplotlib.colors import LogNorm

doFITS=myFITS()
doHVC=seachhvc()


#dataPath="/share/data/mwisp/R19/"

#dataPath="/home/qzyan/WORK/projects/searchHVC/data/G210/"

dataPath="/share/data/mwisp/G130+00/"

if 0:

    allCOFITS=doHVC.getAllCO12FITS(dataPath)

    for eachCO in allCOFITS:
        print eachCO
        doHVC.searchCloud(eachCO,outPath=doHVC.tmpPath )

if 0:
    comBinTB=doHVC.getAllTBFiles(doHVC.tmpPath )

    #print len( comBinTB  )

if 0:
    doHVC.checkTB("combinedTB.fit")

if 1:

    TB=Table.read("combinedTB.fit")


    part1TB= doFITS.selectTBByColRange(TB,"v_cen", minV=30 )

    doHVC.drawCloudSpectra(part1TB ,doHVC.tmpPath )


    part2TB= doFITS.selectTBByColRange(TB,"v_cen", maxV=-90 )
    doHVC.drawCloudSpectra(part2TB ,doHVC.tmpPath)
