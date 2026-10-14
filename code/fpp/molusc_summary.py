# Copied from ../../211038138/fpp_joint/scripts/molusc_summary.py with M=0.381 (TIC 118598568), plx 6.7936 mas, paths for EPIC 247455141.
import pandas as pd, numpy as np, glob
F='k2-deep/247455141/fpp_joint/molusc_out/'
fs=sorted(glob.glob(F+'m100k_s2?_all.csv')); A=pd.concat([pd.read_csv(f) for f in fs],ignore_index=True)
n=len(A); M=0.381; A['m2']=A['mass ratio']*M; kept=A[A['Full Rejected']==0]
print('files',len(fs),'generated',n,'kept',len(kept),'frac %.4f'%(len(kept)/n))
print('reject: Gaia %.4f RUWE %.4f'%(A['Gaia Rejected'].mean(), A['RUWE Rejected'].mean()))
st=A['m2']>=0.08
print('stellar (m2>=0.08) generated',st.sum(),'kept',(st&(A['Full Rejected']==0)).sum(),'frac %.4f'%((st&(A['Full Rejected']==0)).sum()/st.sum()))
for f in fs:
    a=pd.read_csv(f); print(f.split('/')[-1],'kept frac %.4f'%(1-a['Full Rejected'].mean()))
sb=[0,1.1,5,20,60,130,250,1e5]; gb=[-0.2,1,2,3,4,5,7.1,99]
B=A[st].copy(); B['s']=pd.cut(B['Projected Separation(AU)'],sb); B['g']=pd.cut(B['DeltaG'],gb)
print('Fraction of stellar companions REJECTED (rows projected sep AU; 1 arcsec = %.1f AU)'%(1000/6.7936))
print(B.groupby(['s','g'],observed=True)['Full Rejected'].mean().unstack().round(2).to_string())
rng=np.random.default_rng(7); N=1000000; nk=int(round(len(kept)/n*N))
kept[['mass ratio','semi-major axis(AU)','eccentricity']].iloc[rng.integers(0,len(kept),nk)].reset_index(drop=True).to_csv(F+'molusc_for_trc_N1e6.csv')
print('wrote',nk,'rows')
