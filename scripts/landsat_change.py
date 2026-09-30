import numpy as np, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pyproj import Transformer
L=np.load('data/landsat.npz'); yrs=list(L['years']); S=L['stack']; T=L['T']; n,_,H,W=S.shape
r,c=np.mgrid[0:H,0:W]; lon,lat=Transformer.from_crs(32641,4326,always_xy=True).transform(T[2]+(c+.5)*30,T[5]-(r+.5)*30)
def sw(k): return (S[k][4]+S[k][5])/2
def vis(k): return (S[k][0]+S[k][1]+S[k][2])/3
res={}
base=yrs.index(2020) if 2020 in yrs else 0
for k,y in enumerate(yrs):
    # change vs 2020 in SWIR brightness, normalised by each scene's median (sensor/illumination robust)
    a=sw(k)/np.nanmedian(sw(k)); b=sw(base)/np.nanmedian(sw(base)); d=a-b
    res[int(y)]=d
json.dump({int(y):None for y in yrs},open('landsat_years.json','w'))
np.savez_compressed('data/landsat_change.npz',**{str(k):v.astype('float32') for k,v in res.items()},lon=lon.astype('float32'),lat=lat.astype('float32'))
def box(a,b,cc,d): return (lon>=a)&(lon<=cc)&(lat>=b)&(lat<=d)
zones={'Given point ±1.5 km':box(62.417,28.958,62.448,28.986),'Western complex':box(62.37,28.97,62.41,29.00),'Control plain S':box(62.30,28.88,62.40,28.92)}
for y in yrs:
    print(y,{n:round(float(np.nanmean(np.abs(res[int(y)][m])>0.15)*100),2) for n,m in zones.items()})
