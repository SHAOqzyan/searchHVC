
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import scipy.odr.odrpack as odrpack
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

from mwispDBSCAN import MWISPDBSCAN
import glob
import matplotlib.pyplot as plt
from myPYTHON import *
from astropy.io import fits
import glob
import os
import sys
import seaborn as sns
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy import optimize
from progressbar import *
from astropy.table import Table,vstack
import gc
import scipy

doFITS=myFITS()
doMWdbscan= MWISPDBSCAN()

class seachhvc(object):

    dataPath="./data/"

    tmpPath = "./tmpFiles/"
    figurePath= "./tmpFigures/"
    COtags=["U","L","L2"]

    #cutVrange=[-1300  ,1110]
    cutVrange=[-600  ,600]

    lRad=0.25
    bRange=0.25
    tbNameCol="tbNameCol"
    def __init__(self):
        pass




    def getLBfromFileName(self,nameStr):

        """
        return
        :param nameStr:
        :return:
        """
        nameStr=os.path.basename(nameStr )
        for eachTag in self.COtags:

            if eachTag in nameStr:
                nameStr = nameStr.split( eachTag )[0]

        splitCor="+"

        if "-" in nameStr:
            splitCor="-"



        lStr,bStr=nameStr.split(splitCor )


        l,b=map(float,[lStr,splitCor+bStr])

        l=l/10.
        b=b/10.

        return l,b


    def getRMSfitsName(self,rawCOfits):

        """
        return the corresponding rms fits
        :param rawCOfits:
        :return:
        """

        dirPath= os.path.dirname(rawCOfits)
        baseName=os.path.basename(rawCOfits )

        name,ext= os.path.splitext( baseName  )

        return os.path.join(dirPath, name+"_rms"+ext )



    def searchCloud(self,rawCOFITS,path=None ):
        """
        run dbscan, down to 2 sigma
        :param rawCOFITS:
        :return:
        """

        if path is not None:

            rawCOFITS=os.path.join(path,rawCOFITS)

        #cropFITS

        centerL,centerB= self.getLBfromFileName(rawCOFITS)

        cutLrange= [centerL-self.lRad, centerL+self.lRad ]

        cutBrange= [centerB-self.lRad, centerB+self.lRad ]

        rmsFITS=self.getRMSfitsName(rawCOFITS)

        rawCOFITS=doFITS.cropFITS(rawCOFITS,Vrange=self.cutVrange,Lrange= cutLrange ,Brange=cutBrange, velUnit ='ms',overWrite=True)
        rmsFITS=doFITS.cropFITS2D(rmsFITS, Lrange= cutLrange ,Brange=cutBrange ,overWrite=True)


        doMWdbscan.rawCOFITS= rawCOFITS
        doMWdbscan.rmsFITS =  rmsFITS
        doMWdbscan.processPath = self.tmpPath

        doMWdbscan.computeDBSCAN()

        doMWdbscan.getCatFromLabelArray(doClean=True)

        doMWdbscan.produceCleanFITS()


    def getAllTBFiles(self,tbPath):
        """

        :param tbPath:
        :return:
        """
        searchStr= os.path.join(tbPath,"*_Clean.fit")



        searchResesults= glob.glob(searchStr )
        #return searchResesults
        combineTB=None


        for eachTB in searchResesults:

            testTB=Table.read(eachTB)

            #add a column of FITSnames

            tbName=os.path.basename(eachTB)

            if len(testTB)==0:
                continue


            testTB[self.tbNameCol]= tbName

            if combineTB is None:
                    combineTB=testTB

            else:
                combineTB=vstack([combineTB,testTB])
        combineTB.write("combinedTB.fit",overwrite=True )
        return combineTB







    def getAllCO12FITS(self,dataPath):

        """
        return all files of  CO12
        :param dataPath:
        :return:
        """

        searchStr= os.path.join(dataPath,"*U.fits")



        searchResesults= glob.glob(searchStr )
        return searchResesults


    def checkTB(self,TBName):

        """
        draw figure of vDistribution

        :param TBName:
        :return:
        """
        TB=Table.read(TBName)

        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        ax.scatter(TB["v_cen"],TB["peak"],s=5,color='blue')


        plt.savefig("checkTB.png", bbox_inches='tight')


    def drawCloudSpectra(self,TB,coPath):
        """
        draw all spectra in a TB
        :param TB:
        :return:
        """
        fig = plt.figure(figsize=(12, 12))
        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        for eachRow in TB:

            rawCOName= eachRow[self.tbNameCol]

            rawCOFITS = os.path.join(  coPath , rawCOName[0:9]+"_C.fits" )

            drawL,drawB = eachRow["x_cen"] , eachRow["y_cen"]

            Radius=30./3600.*1


            lRangeSpec = [ drawL - Radius, drawL + Radius  ]
            bRangeSpec = [ drawB - Radius, drawB + Radius  ]

            averageSpec,Vs= doFITS.getAverageSpecByLBrange(rawCOFITS,  lRange = lRangeSpec  ,bRange= bRangeSpec )

            spectralMean=np.mean(  averageSpec[ int(eachRow["peakV"])-50:int(eachRow["peakV"]) -10  ]  )

            rms=np.sqrt( np.mean(  np.square(averageSpec[averageSpec<= spectralMean ]) )  )
            #ax.axhline(Vs, averageSpec,color='blue',where='mid',lw=0.8   ,zorder=2, label= labelStr )

            ######get goo channels

            cutSpec= averageSpec[ int(eachRow["peakV"])-10:int(eachRow["peakV"])+10  ]

            SNR3Pix=cutSpec[cutSpec>= spectralMean+rms*3 ]


            if len( SNR3Pix )<=2 :
                #insufficient dispersion


                continue

            ax = fig.add_subplot(2, 1, 1)
            ax2 = fig.add_subplot(2, 1, 2, sharex=ax)

            ax.axhline(y=spectralMean+rms*3, ls="--", color='black', lw=0.8)

            ax.plot([eachRow["v_cen"], eachRow["v_cen"]], [-0.5, np.nanmax(averageSpec) * 1.2], color='red', lw=1.5,
                    zorder=1)
            ax2.plot([eachRow["v_cen"], eachRow["v_cen"]], [-0.5, 5], color='red', lw=1.5, zorder=1)

            peakL, peakB ,peakV= eachRow["peakL"], eachRow["peakB"], eachRow["peakV"]

            labelStr = "peak index(L,B,V): ({}, {}, {}  )".format(int(peakL), int(peakB), int(peakV) )

            ax.step(Vs, averageSpec, color='blue', where='mid', lw=0.8, zorder=2, label=labelStr)


            ax.axhline(y= spectralMean , ls="-", color='black', lw=0.8)
            ax2.scatter(TB["v_cen"],  TB["peak"],color='green' ,s=5 )
            ax.set_xlim(eachRow["v_cen"]-20, eachRow["v_cen"]+20   )


            ax.legend(loc=1,handlelength=0.5)
            saveFigname= rawCOName[0:9]+"{}.png".format( eachRow["_idx"] )

            saveFigname= os.path.join( self.figurePath, saveFigname  )

            ax.set_xlabel("Radial velocity")
            ax.set_ylabel(r"Average T$_{\rm mb}$ (5$\times$5)")
            ax2.set_ylabel(r"Peak T$_{\rm mb}$")


            plt.savefig( saveFigname , bbox_inches='tight',dpi=100)

            ax.cla()
            ax2.cla()




    def ZZZ(self):
        pass