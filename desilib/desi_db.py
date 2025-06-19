import fitsio
import sys
import os
import pandas as pd
import numpy
#pylibdir=os.environ['PYLIB']
#sys.path.append(pylibdir)
#import sdss_catalog

#spall=sdss_catalog.spall()
#print(spall.dr)
#print(spall.fitstablename)
#spall.read()

class zall():
      def __init__(self):
          self.version='iron'
          self.dr='DR1'
          self.fitstablename=os.environ['DESI_REDUX']+'zall-tilecumulative-'+self.version+'.fits'
          h=fitsio.read_header(self.fitstablename,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)

      def readstar(self):
          columns=['PLATE','MJD','FIBERID','OBJTYPE','CLASS'=='STAR','SUBCLASS']
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
          self.platelist=d['PLATE']
          self.mjdlist=d['MJD']
          self.fiberlist=d['FIBERID']
          self.objtypelist=d['OBJTYPE']
          self.classlist=d['CLASS']
          self.subclasslist=d['SUBCLASS']

      def readtest(self):
          columns=['TARGETID','SURVEY','PROGRAM',\
          'FIRSTNIGHT','LASTNIGHT',\
          'Z','ZERR','ZWARN',\
          'TILEID','TARGET_RA','TARGET_DEC','SPECTYPE','PRIORITY',\
          'DEVICE_LOC','PETAL_LOC','LOCATION','FIBER']
          print('Reading ',self.fitstablename)
#          print('Number of Spectra is ',len(self.nspec))
          d=fitsio.read(self.fitstablename,columns=columns,rows=self.rows)
#          d=fitsio.read(self.fitstablename,columns=columns,rows=1000)
          self.targetid=d['TARGETID']
          self.survey=d['SURVEY']
          self.program=d['PROGRAM']

          self.firstnight=d['FIRSTNIGHT']
          self.lastnight=d['LASTNIGHT']

          self.z=d['Z']
          self.zerr=d['ZERR']
          self.zwarn=d['ZWARN']

          self.tileid=d['TILEID']
          self.ra=d['TARGET_RA']
          self.dec=d['TARGET_DEC']
          self.spectype=d['SPECTYPE']
          self.device_loc=d['DEVICE_LOC']
          self.petal_loc=d['PETAL_LOC']
          self.priority=d['PRIORITY']
          self.location=d['LOCATION']
          self.fiber=d['FIBER']
          #self.subclass=d['SUBCLASS']
          #self.specprimary=d['SPECPRIMARY']

zall=zall()
print(zall.nspec)
zall.readtest()
print('targetid',zall.targetid)
print('survey',zall.survey)
print('program',zall.program) 

print('firstnight',zall.firstnight) 
print('lastnight',zall.lastnight) 

print('z',zall.z)
print('zerr',zall.zerr)
print('zwarn',zall.zwarn)

print('tileid',zall.tileid)
print('ra',zall.ra)
print('dec',zall.dec)
print('priority',zall.priority)
print('spectype',zall.spectype)
print('location',zall.location)
print('device_loc',zall.device_loc)
print('petal_loc',zall.petal_loc)
print('fiber',zall.fiber)
#print('tileid',zall.tileid)
