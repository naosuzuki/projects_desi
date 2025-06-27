import pandas as pd
import os

def extract_stars(csvfile):
   df=pd.read_csv(csvfile,dtype={'targetid':int,'ra':'float32','dec':'float32','spectype':str})
   print(len(df))

   dfstar1=df[(df['spectype']=='STAR') & (df['targetid']>0)]
   print(len(dfstar1))
#   dfstar2=dfstar1[['targetid','ra','dec']]
#   dfstar2[0:2000000].to_csv('desistar1_float32.csv',index=False)
#   dfstar2[2000001:].to_csv('desistar2_float32.csv',index=False)

# DR2
csvfile=os.environ['DESI_REDUX']+'zall-tilecumulative-loa.csv'
# DR1
csvfile=os.environ['DESI_REDUX']+'zall-tilecumulative-iron.csv'
print(csvfile)
extract_stars(csvfile)
