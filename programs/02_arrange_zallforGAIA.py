import pandas as pd
import os

def extract_stars(csvfile,output):
   df=pd.read_csv(csvfile,dtype={'targetid':int,'ra':'float32','dec':'float32','spectype':str,'objtype':str})
   print(len(df))

   dfstar1=df[(df['spectype']=='STAR') & (df['targetid']>0) & (df['objtype']=='TGT')]
   print(len(dfstar1))
   dfstar2=dfstar1[['targetid','ra','dec']]
   #dfstar2.to_csv(output,index=False)
   dfstar2[0:3000000].to_csv('desidr2pixstar1_float32.csv',index=False)
   dfstar2[3000001:6000000].to_csv('desidr2pixstar2_float32.csv',index=False)
   dfstar2[6000001:9000000].to_csv('desidr2pixstar3_float32.csv',index=False)
   dfstar2[9000001:12000000].to_csv('desidr2pixstar4_float32.csv',index=False)
   dfstar2[12000001:].to_csv('desidr2pixstar5_float32.csv',index=False)

# DR1
csvfile=os.environ['DESI_REDUX']+'zall-tilecumulative-iron.csv'
output='../csvfiles/desidr1star_float32.csv'
# DR2
csvfile=os.environ['DESI_REDUX']+'zall-tilecumulative-loa.csv'
output='../csvfiles/desidr2star_float32.csv'
csvfile=os.environ['DESI_REDUX']+'zall-pix-loa.csv'
output='../csvfiles/desidr2pixstar_float32.csv'
print(csvfile)
extract_stars(csvfile,output)
