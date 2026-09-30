import numpy as np, rasterio, glob, re
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
z=np.load('data/emit_feat.npz'); glon,glat=z['glon'],z['glat']; H,W=len(glat),len(glon)
T=from_origin(glon[0]-0.00025,glat[0]+0.00025,0.0005,0.0005)
grans=sorted({re.search(r'AST_07XT_(\d+)_',f).group(1) for f in glob.glob('data/aster/AST_07XT_*')})
acc={k:[] for k in ('AlOH','kaol','ferric')}
for g in grans:
    B={}; bad=[]
    for b in ['01','02','03N','04','05','06','07','08','09']:
        fs=glob.glob(f'data/aster/AST_07XT_{g}_*B{b}.tif')
        if not fs: bad.append(b); continue
        a=np.full((H,W),np.nan,'float32')
        with rasterio.open(fs[0]) as s:
            src=s.read(1).astype('float32'); nod=src<=0; src=src*0.001; src[nod]=np.nan; v=src[np.isfinite(src)]
            if v.size and np.median(v)>=0.999*v.max(): bad.append(b+'(clipped)')
            reproject(src,a,src_transform=s.transform,src_crs=s.crs,dst_transform=T,dst_crs='EPSG:4326',resampling=Resampling.average,src_nodata=np.nan,dst_nodata=np.nan)
        B[b]=a
    if bad: print(g,'REJECT',bad); continue
    ndvi=(B['03N']-B['02'])/(B['03N']+B['02'])
    ok=np.all([np.isfinite(B[b])&(B[b]>0.02) for b in B],0)&(ndvi<0.25)
    print(g,'usable',round(float(ok.mean()),3))
    acc['AlOH'].append(np.where(ok,(B['05']+B['07'])/B['06'],np.nan)); acc['kaol'].append(np.where(ok,(B['04']+B['06'])/B['05'],np.nan)); acc['ferric'].append(np.where(ok,B['02']/B['01'],np.nan))
out={k:np.nanmedian(np.stack(v),0).astype('float32') for k,v in acc.items()}
np.savez_compressed('data/aster_idx.npz',**out)
def rz(a): med=np.nanmedian(a); mad=1.4826*np.nanmedian(abs(a-med)); return (a-med)/mad
F=np.load('data/facies.npz'); LON,LAT=np.meshgrid(glon,glat); ctrl=(LON>=62.30)&(LON<=62.40)&(LAT>=28.88)&(LAT<=28.92)
A=rz(out['AlOH']); K=rz(out['kaol'])
for n,m in (('EMIT phyllic zones',F['phy']),('EMIT argillic zones',F['arg']),('EMIT leached cap',F['cap']),('Alluvial control',ctrl)):
    print(f'{n:22s} ASTER AlOH z {np.nanmedian(A[m]):+.2f}   ASTER kaol z {np.nanmedian(K[m]):+.2f}   n={int(m.sum())}')
e=np.nanmean([rz(z[f'{d}_AlOH_2200']) for d in ('20231001','20241001')],0); v=np.isfinite(e)&np.isfinite(A)
print('pixelwise r(EMIT AlOH, ASTER AlOH) =',round(float(np.corrcoef(e[v],A[v])[0,1]),3))
