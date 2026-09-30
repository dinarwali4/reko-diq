import numpy as np, json, rasterio
from rasterio.merge import merge
from scipy import ndimage as ndi
from pyproj import Geod
z=np.load('data/emit_feat.npz'); glon,glat=z['glon'],z['glat']; LON,LAT=np.meshgrid(glon,glat)
def rz(a): med=np.nanmedian(a); mad=1.4826*np.nanmedian(abs(a-med)); return (a-med)/mad
dates=['20231001','20241001']
hits=[];S=[];K=[];A=[];J=[];F=[]
for d in dates:
    a,k,j,f=rz(z[f'{d}_AlOH_2200']),rz(z[f'{d}_Kaol_2165']),rz(z[f'{d}_Jar_2265']),rz(z[f'{d}_Fe3_900'])
    n=(np.nan_to_num(a)>2).astype(int)+(np.nan_to_num(k)>2)+(np.nan_to_num(j)>2)
    hits.append(n>=2); S.append((np.nan_to_num(a)+np.nan_to_num(k)+np.nan_to_num(j))/3); A.append(a);K.append(k);J.append(j);F.append(f)
zone=ndi.binary_opening(hits[0]&hits[1],structure=np.ones((2,2)))
sc=np.mean(S,0); Am,Km,Jm,Fm=[np.nanmean(x,0) for x in (A,K,J,F)]
src=[rasterio.open(f'data/dem_N{t}E062.tif') for t in (28,29)]; mos,_=merge(src,bounds=(62.30,28.88,62.55,29.06),res=(0.0005,0.0005)); D=mos[0][:len(glat),:len(glon)]
gy,gx=np.gradient(D.astype(float),55.3,48.6); slope=np.degrees(np.arctan(np.hypot(gx,gy)))
# control: alluvial plains south (no outcrop) 62.30-62.40, 28.88-28.92
ctrl=(LON>=62.30)&(LON<=62.40)&(LAT>=28.88)&(LAT<=28.92)
print('control cells in zone:',int((zone&ctrl).sum()),'of',int(ctrl.sum()))
lab,n=ndi.label(ndi.binary_dilation(zone,iterations=1)&zone,structure=np.ones((3,3))); sz=ndi.sum(zone,lab,range(1,n+1))
G=Geod(ellps='WGS84'); P=(28.972036,62.432438); pix_ha=0.0005*111320*np.cos(np.radians(28.97))*0.0005*110574/1e4
Zs=[]
for kk in np.argsort(-sz)+1:
    m=lab==kk
    if m.sum()<4: continue
    w=np.where(m,sc,-9); i=np.unravel_index(np.argmax(w),w.shape); _,_,dist=G.inv(P[1],P[0],LON[i],LAT[i])
    dom=max((('phyllic (sericite)',np.nanmean(Am[m])),('argillic (kaolinite/alunite)',np.nanmean(Km[m])),('leached cap (jarosite)',np.nanmean(Jm[m]))),key=lambda t:t[1])[0]
    Zs.append(dict(cells=int(m.sum()),area_ha=round(float(m.sum()*pix_ha),1),peak=[round(float(LAT[i]),5),round(float(LON[i]),5)],centroid=[round(float(LAT[m].mean()),5),round(float(LON[m].mean()),5)],
      score=round(float(sc[i]),2),z_sericite=round(float(np.nanmean(Am[m])),2),z_kaol=round(float(np.nanmean(Km[m])),2),z_jarosite=round(float(np.nanmean(Jm[m])),2),z_ferric=round(float(np.nanmean(Fm[m])),2),
      dominant=dom,elev=[int(D[m].min()),int(D[m].max())],slope=round(float(slope[m].mean()),1),dist_given_km=round(dist/1000,2)))
for j,t in enumerate(Zs): t['id']=f'RZ{j+1}'
for t in Zs[:18]: print(t)
json.dump(Zs,open('zones.json','w'),indent=1)
np.savez_compressed('data/zone_grid.npz',zone=zone,score=sc.astype('float32'),Am=Am.astype('float32'),Km=Km.astype('float32'),Jm=Jm.astype('float32'),Fm=Fm.astype('float32'),glon=glon,glat=glat,D=D)
