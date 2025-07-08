import pandas as pd
import os
import sys

desilibdir=os.environ['DESILIB']
sys.path.append(desilibdir)
import desi_db

zall=desi_db.zall()
zall.read_pix()

df=pd.DataFrame({'targetid':zall.targetid,\
                 'ra':zall.target_ra,\
                 'dec':zall.target_dec,\
                 'spectype':zall.spectype,\
                 'healpix':zall.healpix,\
                 'survey':zall.survey,\
                 'program':zall.program,\
                 'spgrval':zall.spgrval,\
                 'coadd_exptime':zall.coadd_exptime})

#spall=sdss_catalog.spall()
#print(spall.dr)
#print(spall.fitstablename)
#spall.read()


