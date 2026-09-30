import numpy as np, json
from scipy import ndimage as ndi
from pyproj import Geod
z=np.load('data/emit_feat.npz'); glon,glat=z['glon'],z['glat']; LON,LAT=np.meshgrid(glon,glat)
G=np.load('data/zone_grid.npz'); D=G['D']
def rz(a): med=np.nanmedian(a); mad=1.4826*np.nanmedian(abs(a-med)); return (a-med)/mad
dates=['20231001','20241001']
def both(fn): return np.logical_and(*[fn(d) for d in dates])
phy=both(lambda d: rz(z[f'{d}_AlOH_2200'])>2.5)
arg=both(lambda d: rz(z[f'{d}_Kaol_2165'])>2.5)
cap=both(lambda d: (rz(z[f'{d}_Jar_2265'])>2.5)&(rz(z[f'{d}_Fe3_900'])>0))
op=lambda m: ndi.binary_opening(m,structure=np.ones((2,2)))
phy,arg,cap=op(phy),op(arg),op(cap)
ctrl=(LON>=62.30)&(LON<=62.40)&(LAT>=28.88)&(LAT<=28.92)
for n,m in (('phyllic',phy),('argillic',arg),('leached cap',cap)): print(n,'cells',int(m.sum()),'control',int((m&ctrl).sum()))
Gd=Geod(ellps='WGS84'); P=(28.972036,62.432438); pix=0.0005*111320*np.cos(np.radians(28.97))*0.0005*110574/1e4
Zs=[]
for name,m in (('phyllic (sericite)',phy),('argillic (kaolinite/alunite)',arg),('leached cap (jarosite + Fe)',cap)):
    lab,n=ndi.label(ndi.binary_dilation(m,iterations=2),structure=np.ones((3,3)))
    for k in range(1,n+1):
        mm=(lab==k)&m
        if mm.sum()<6: continue
        la,lo=float(LAT[mm].mean()),float(LON[mm].mean()); _,_,d=Gd.inv(P[1],P[0],lo,la)
        Zs.append(dict(facies=name,cells=int(mm.sum()),area_ha=round(float(mm.sum()*pix),1),centroid=[round(la,5),round(lo,5)],dist_given_km=round(d/1000,2),elev=[int(D[mm].min()),int(D[mm].max())]))
Zs=sorted(Zs,key=lambda t:-t['area_ha'])
for j,t in enumerate(Zs): t['id']=f'RZ{j+1}'; print(t)
json.dump(Zs,open('facies_zones.json','w'),indent=1)
np.savez_compressed('data/facies.npz',phy=phy,arg=arg,cap=cap,glon=glon,glat=glat)
