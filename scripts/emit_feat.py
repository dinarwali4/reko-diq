import numpy as np, glob
feat={'Kaol_2165':([2165],2120,2185),'Fe3_900':([900],760,1250),'AlOH_2200':([2195,2205],2140,2255),'Jar_2265':([2262],2240,2290),'MgOH_2330':([2325,2335],2260,2390)}
out={}
for a in sorted(glob.glob('data/emit_aoi_*.npz')):
    d=a.split('_')[-1][:8]; z=np.load(a); C=z['cube'].astype('float32'); wl=z['wl']; mk=z['mask']
    bi=lambda w:int(np.argmin(abs(wl-w)))
    r671,r800,r560,r1610=[C[...,bi(w)] for w in (671,800,560,1610)]
    ndvi=(r800-r671)/(r800+r671); ndsi=(r560-r1610)/(r560+r1610)
    valid=(mk==0)&np.isfinite(r671)&(r671>0.03)&(ndvi<0.30)&(ndsi<0.35)
    res={'valid':valid,'ndvi':ndvi.astype('float32'),'ndsi':ndsi.astype('float32')}
    for n,(c,l,r) in feat.items():
        Rc=np.nanmean([C[...,bi(w)] for w in c],0); L=C[...,bi(l)]; Rr=C[...,bi(r)]
        t=(np.mean([wl[bi(w)] for w in c])-wl[bi(l)])/(wl[bi(r)]-wl[bi(l)]); dd=1-Rc/(L+(Rr-L)*t); dd[~valid]=np.nan; res[n]=dd.astype('float32')
    out[d]=res; print(d,'valid',round(float(valid.mean()),3),'snow',round(float(np.nanmean(ndsi>0.35)),3),{n:round(float(np.nanmedian(res[n])),4) for n in feat})
np.savez_compressed('data/emit_feat.npz',dates=np.array(list(out)),glon=z['glon'],glat=z['glat'],**{f'{d}_{k}':v for d,r in out.items() for k,v in r.items()})
