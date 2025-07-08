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
                 'spgrpval':zall.spgrpval,\
		 'coadd_numexp':zall.coadd_numexp,\
		 'coadd_numnight':zall.coadd_numnight,\
		 'coadd_numtile':zall.coadd_numtile,\
                 'coadd_exptime':zall.coadd_exptime})

df.to_csv('zalldr2pix.csv',index=False)

#spall=sdss_catalog.spall()
#print(spall.dr)
#print(spall.fitstablename)
#spall.read()


