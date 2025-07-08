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
                 'objtype':zall.objtype,\
                 'gaia_g':zall.gaia_g,\
		 'coadd_numexp':zall.coadd_numexp,\
		 'coadd_numnight':zall.coadd_numnight,\
		 'coadd_numtile':zall.coadd_numtile,\
                 'coadd_exptime':zall.coadd_exptime})

df.to_csv('zalldr2pix.csv',index=False)

del df

df=pd.read_csv('zalldr2pix.csv')
dfstar=df[df['spectype']=='STAR']
dfstar.to_csv('dfstar.csv',index=False)

sys.exit(1)

#df['spectype'] = df['spectype'].astype(str)
#dfstar=df[df['spectype']=='STAR']
#for col in df.columns:
#    print(col, df[col].dtype)

# Force the 'spectype' column to be native-endian string
#df['spectype'] = df['spectype'].astype(str)
#df['survey'] = df['survey'].astype(str)
#df['program'] = df['program'].astype(str)
#df['objtype'] = df['objtype'].astype(str)

#for col in df.columns:
#    print(col, df[col].dtype)
# Now filtering works safely
#dfstar=df[df['spectype']=='STAR']
#dfstar.to_csv('zalldr2pixstar.csv',index=False)

#spall=sdss_catalog.spall()
#print(spall.dr)
#print(spall.fitstablename)
#spall.read()


