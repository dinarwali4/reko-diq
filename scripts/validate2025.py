import numpy as np, json
z=np.load('data/emit_feat.npz'); glon,glat=z['glon'],z['glat']; LON,LAT=np.meshgrid(glon,glat); d='20250202'
def rz(a): med=np.nanmedian(a); mad=1.4826*np.nanmedian(abs(a-med)); return (a-med)/mad
F=np.load('data/facies.npz'); ctrl=(LON>=62.30)&(LON<=62.40)&(LAT>=28.88)&(LAT<=28.92)
A,K,J,Fe=[rz(z[f'{d}_{k}']) for k in ('AlOH_2200','Kaol_2165','Jar_2265','Fe3_900')]
res={}
for n,m,key,arr,thr in (('Phyllic zones',F['phy'],'sericite z>2.5',A,2.5),('Argillic zones',F['arg'],'kaolinite z>2.5',K,2.5),('Leached-cap zones',F['cap'],'jarosite z>2.5',J,2.5)):
    v=np.isfinite(arr)
    res[n]=dict(frac_in_zone=round(float(np.nanmean(arr[m&v]>thr)*100),1),frac_in_control=round(float(np.nanmean(arr[ctrl&v]>thr)*100),2),median_z_zone=round(float(np.nanmedian(arr[m&v])),2),median_z_control=round(float(np.nanmedian(arr[ctrl&v])),2),valid=round(float(v[m].mean()),2))
    print(n,res[n])
json.dump(res,open('validation_2025.json','w'),indent=1)
