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
          self.fitstablename_tile=os.environ['DESI_REDUX']+'zall-tilecumulative-'+self.version+'.fits'
          self.fitstablename_pix=os.environ['DESI_REDUX']+'zall-pix-'+self.version+'.fits'

      def read_pix(self):
          columns=['TARGETID','SURVEY','PROGRAM','HEALPIX','SPGRPVAL',\
          'Z','ZERR','ZWARN','CHI2',\
          'SPECTYPE','SUBTYPE',\
          'TARGET_RA','TARGET_DEC','OBJTYPE',\
          'FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2',\
          'FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2',\
          'GAIA_PHOT_G_MEAN_MAG','GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG',\
          'PHOTSYS','PARALLAX',\
          'COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT','COADD_NUMTILE']
          h=fitsio.read_header(self.fitstablename_pix,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)
          d=fitsio.read(self.fitstablename_pix,columns=columns,rows=self.rows)
          self.targetid=d['TARGETID']
          self.lastnight=d['LASTNIGHT']
          self.z=d['Z']
          self.zerr=d['ZERR']
          self.zwarn=d['ZWARN']
          self.spectype=d['SPECTYPE']
          self.subtype=d['SUBTYPE']
          self.tileid=d['TILEID']
          self.petal_loc=d['PETAL_LOC']
          self.fiber=d['FIBER']
          self.target_ra=d['TARGET_RA']
          self.target_dec=d['TARGET_DEC']
          self.fiber_ra=d['MEAN_FIBER_RA']
          self.fiber_dec=d['MEAN_FIBER_DEC']
          self.fiberstd_ra=d['STD_FIBER_RA']
          self.fiberstd_dec=d['STD_FIBER_DEC']
          self.photsys=d['PHOTSYS']
          self.parallax=d['PARALLAX']
          self.flux_g=d['FLUX_G']
          self.flux_r=d['FLUX_R']
          self.flux_z=d['FLUX_Z']
          self.flux_w1=d['FLUX_W1']
          self.flux_w2=d['FLUX_W2']
          self.flux_ivar_g=d['FLUX_IVAR_G']
          self.flux_ivar_r=d['FLUX_IVAR_R']
          self.flux_ivar_z=d['FLUX_IVAR_Z']
          self.flux_ivar_w1=d['FLUX_IVAR_W1']
          self.flux_ivar_w2=d['FLUX_IVAR_W2']
          self.gaia_g=d['GAIA_PHOT_G_MEAN_MAG']
          self.gaia_bp=d['GAIA_PHOT_BP_MEAN_MAG']
          self.gaia_rp=d['GAIA_PHOT_RP_MEAN_MAG']
          self.coadd_numexp=d['COADD_NUMEXP']
          self.coadd_exptime=d['COADD_EXPTIME']

      def read_tile(self):
          columns=['TARGETID','LASTNIGHT','Z','ZERR','ZWARN','CHI2',\
          'SPECTYPE','SUBTYPE',\
          'TILEID','PETAL_LOC','FIBER','TARGET_RA','TARGET_DEC',\
          'MEAN_FIBER_RA','STD_FIBER_RA',\
          'MEAN_FIBER_DEC','STD_FIBER_DEC',\
          'PHOTSYS','PARALLAX',\
          'FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2',\
          'FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2',\
          'GAIA_PHOT_G_MEAN_MAG','GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG',\
          'COADD_NUMEXP','COADD_EXPTIME']
          h=fitsio.read_header(self.fitstablename_tile,ext=1)
          self.nspec=h['NAXIS2']
          self.rows=numpy.arange(self.nspec)
          d=fitsio.read(self.fitstablename_tile,columns=columns,rows=self.rows)
          self.targetid=d['TARGETID']
          self.lastnight=d['LASTNIGHT']
          self.z=d['Z']
          self.zerr=d['ZERR']
          self.zwarn=d['ZWARN']
          self.spectype=d['SPECTYPE']
          self.subtype=d['SUBTYPE']
          self.tileid=d['TILEID']
          self.petal_loc=d['PETAL_LOC']
          self.fiber=d['FIBER']
          self.target_ra=d['TARGET_RA']
          self.target_dec=d['TARGET_DEC']
          self.fiber_ra=d['MEAN_FIBER_RA']
          self.fiber_dec=d['MEAN_FIBER_DEC']
          self.fiberstd_ra=d['STD_FIBER_RA']
          self.fiberstd_dec=d['STD_FIBER_DEC']
          self.photsys=d['PHOTSYS']
          self.parallax=d['PARALLAX']
          self.flux_g=d['FLUX_G']
          self.flux_r=d['FLUX_R']
          self.flux_z=d['FLUX_Z']
          self.flux_w1=d['FLUX_W1']
          self.flux_w2=d['FLUX_W2']
          self.flux_ivar_g=d['FLUX_IVAR_G']
          self.flux_ivar_r=d['FLUX_IVAR_R']
          self.flux_ivar_z=d['FLUX_IVAR_Z']
          self.flux_ivar_w1=d['FLUX_IVAR_W1']
          self.flux_ivar_w2=d['FLUX_IVAR_W2']
          self.gaia_g=d['GAIA_PHOT_G_MEAN_MAG']
          self.gaia_bp=d['GAIA_PHOT_BP_MEAN_MAG']
          self.gaia_rp=d['GAIA_PHOT_RP_MEAN_MAG']
          self.coadd_numexp=d['COADD_NUMEXP']
          self.coadd_exptime=d['COADD_EXPTIME']

      def readtest(self):
          columns=['TARGETID','SURVEY','PROGRAM',\
          'FIRSTNIGHT','LASTNIGHT',\
          'Z','ZERR','ZWARN','CHI2','NPIXELS','LAMBDA_REF',\
          'SPECTYPE','SUBTYPE','PRIORITY',\
          'TILEID','TARGET_RA','TARGET_DEC',
          'DEVICE_LOC','PETAL_LOC','LOCATION','FIBER',\
          'REF_ID','REF_CAT',\
          'PHOTSYS','PARALLAX',\
          'GAIA_PHOT_G_MEAN_MAG','GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG',\
          'COADD_NUMEXP','COADD_EXPTIME']
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
          self.chisq=d['CHI2']
          self.npix=d['NPIXELS']

          self.spectype=d['SPECTYPE']
          self.subtype=d['SUBTYPE']

          self.lambda_ref=d['LAMBDA_REF']
          self.tileid=d['TILEID']
          self.ra=d['TARGET_RA']
          self.dec=d['TARGET_DEC']
          self.device_loc=d['DEVICE_LOC']
          self.petal_loc=d['PETAL_LOC']
          self.priority=d['PRIORITY']
          self.location=d['LOCATION']
          self.fiber=d['FIBER']
          self.ref_id=d['REF_ID']
          self.ref_cat=d['REF_CAT']
          self.photsys=d['PHOTSYS']
          self.parallax=d['PARALLAX']
          self.gaia_g=d['GAIA_PHOT_G_MEAN_MAG']
          self.gaia_bp=d['GAIA_PHOT_BP_MEAN_MAG']
          self.gaia_rp=d['GAIA_PHOT_RP_MEAN_MAG']
          self.coadd_numexp=d['COADD_NUMEXP']
          self.coadd_exptime=d['COADD_EXPTIME']

zall=zall()
print(zall.nspec)
zall.readstar()
df=pd.DataFrame({'targetid':zall.targetid,\
                 'ra':zall.target_ra,\
                 'dec':zall.target_dec,\
                 'spectype':zall.spectype,\
                 'tileid':zall.tileid,\
                 'lastnight':zall.lastnight,\
                 'petal_loc':zall.petal_loc,\
                 'fiber':zall.fiber,\
                 'coadd_exptime':zall.coadd_exptime})
df.to_csv('test.csv',index=False)
#print(zall.targetid)
#print(len(zall.targetid))
#print(len(zall.z))
sys.exit(1)
#zall.readtest()
print('gaia_g',zall.gaia_g)
print('gaia_bp',zall.gaia_bp)
print('gaia_rp',zall.gaia_rp)
print('photsys',zall.photsys)
print('parallax',zall.parallax)
print('targetid',zall.targetid)
print('survey',zall.survey)
print('program',zall.program) 

print('firstnight',zall.firstnight) 
print('lastnight',zall.lastnight) 

print('z',zall.z)
print('zerr',zall.zerr)
print('zwarn',zall.zwarn)
print('chi2',zall.chisq)
print('npixels',zall.npix)
print('lambda_ref',zall.lambda_ref)

print('tileid',zall.tileid)
print('ra',zall.ra)
print('dec',zall.dec)
print('priority',zall.priority)
print('spectype',zall.spectype)
print('Unique Object Type',numpy.unique(zall.spectype))
print('subtype',zall.subtype)
print('Unique Subtype',numpy.unique(zall.subtype))
print('location',zall.location)
print('device_loc',zall.device_loc)
print('petal_loc',zall.petal_loc)
print('fiber',zall.fiber)

print('ref_id',zall.ref_id)
print('ref_cat',zall.ref_cat)
print('numexp',zall.coadd_numexp)
print('exptime',zall.coadd_exptime)
#print('tileid',zall.tileid)
