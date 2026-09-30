import numpy as np, json
from scipy import ndimage as ndi
from pyproj import Transformer
f=np.load('data/s2_raw.npz'); S=f['stack'].astype('float32'); ids=list(f['ids']); bn=list(f['bandnames']); B=f['bounds']
dates=[i.split('_')[2] for i in ids]; print(sorted(dates))
scl=S[:,bn.index('SCL')]; good=np.isin(scl,[4,5,7])&(S[:,bn.index('B04')]>0)
R=S[:,:10]/10000.0; R[~np.repeat(good[:,None],10,1)]=np.nan
early=np.array([d<'20250901' for d in dates]); late=np.array([d>='20260601' for d in dates])
print('early',early.sum(),'late',late.sum())
def comp(sel):
    c=np.nanmedian(R[sel],0); b={k:c[i] for i,k in enumerate(bn[:10])}; return b
e,l=comp(early),comp(late)
def bright(b): return (b['B04']+b['B03']+b['B02'])/3
be,bl=bright(e),bright(l); de=be/np.nanmedian(be); dl=bl/np.nanmedian(bl)
d=dl-de
# new bright / graded ground: brightness up >15% relative, persistent in all late scenes
chg=(d>0.15)
chg=ndi.binary_opening(chg,structure=np.ones((3,3)))
lab,n=ndi.label(chg); sz=ndi.sum(chg,lab,range(1,n+1))
H,W=d.shape; tr=Transformer.from_crs(32641,4326,always_xy=True)
out=[]
for k in np.argsort(-sz)[:15]+1:
    m=lab==k; rr,cc=np.where(m); lo,la=tr.transform(B[0]+(cc.mean()+.5)*20,B[3]-(rr.mean()+.5)*20)
    out.append(dict(ha=round(float(m.sum()*0.04),1),centroid=[round(la,5),round(lo,5)],mean_dbright=round(float(d[m].mean()),2)))
for o in out: print(o)
json.dump(out,open('s2_change.json','w'),indent=1)
np.savez_compressed('data/s2_change.npz',d=d.astype('float32'),chg=chg,rgb_late=np.stack([l['B04'],l['B03'],l['B02']]).astype('float32'),rgb_early=np.stack([e['B04'],e['B03'],e['B02']]).astype('float32'),B=B)
