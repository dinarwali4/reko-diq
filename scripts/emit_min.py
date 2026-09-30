import numpy as np, netCDF4, glob, re, json, collections
AOI=(62.30,28.88,62.55,29.06); RES=0.0006
nx=int((AOI[2]-AOI[0])/RES)+1; ny=int((AOI[3]-AOI[1])/RES)+1
fam={'cu_mineral':r'malachite|azurite|chrysocol|cuprite|brochantite|chalcanthite','sulphate':r'jarosit|copiap|schwertmann',
 'ferric':r'hematite|goethite|fe-hydrox|ferrihyd|lepidocro','whitemica':r'muscov|illite|sericit','kaolinite':r'kaolin|dickit|halloys','alunite':r'alunit','pyrophyllite':r'pyrophyl','smectite':r'montmor|beidell|nontron',
 'propylitic':r'chlorit|epidot|prehnite|clinochlore|actinol','carbonate':r'calcite|dolomite|aragonite|siderite|magnesite','nonmineral':r'grass|plastic|tarp|plywood|paint|gas|wood'}
dates=[]; pres={}; hits=[]; stats={}
for f in sorted(glob.glob('data/emit/EMIT_L2B_MIN_001_*.nc')):
    d=f.split('_')[4][:8]; ds=netCDF4.Dataset(f); ds.set_auto_mask(False)
    gt=ds.getncattr('geotransform'); gx=ds['location/glt_x'][:]; gy=ds['location/glt_y'][:]
    names=[str(n) for n in ds['mineral_metadata/name'][:]]; id2n=dict(zip(ds['mineral_metadata/index'][:].tolist(),names))
    H,W=gx.shape; lon=gt[0]+(np.arange(W)+.5)*gt[1]; lat=gt[3]+(np.arange(H)+.5)*gt[5]; LON,LAT=np.meshgrid(lon,lat)
    v=(gx>0)&(gy>0); m=(LON>=AOI[0])&(LON<=AOI[2])&(LAT>=AOI[1])&(LAT<=AOI[3])
    if (m&v).sum()==0: print(d,'no coverage'); continue
    P={k:np.zeros((ny,nx),bool) for k in fam}; retrieved=0; ident=0
    for g in (1,2):
        raw=ds[f'group_{g}_mineral_id'][:]; bd=ds[f'group_{g}_band_depth'][:]
        o=np.full(gx.shape,-9999,'int32'); o[v]=raw[gy[v]-1,gx[v]-1]; b=np.full(gx.shape,np.nan); b[v]=bd[gy[v]-1,gx[v]-1]
        sel=m&v; retrieved=max(retrieved,int(((o!=-9999)&sel).sum())); ident+=int(((o>0)&sel).sum())
        for la,lo,x,bb in zip(LAT[sel&(o>0)],LON[sel&(o>0)],o[sel&(o>0)],b[sel&(o>0)]):
            n=id2n.get(int(x),''); i=int((AOI[3]-la)/RES); j=int((lo-AOI[0])/RES)
            for k,rx in fam.items():
                if re.search(rx,n,re.I): P[k][i,j]=True
            if re.search(fam['cu_mineral']+'|'+fam['sulphate'],n,re.I): hits.append((d,g,round(float(la),5),round(float(lo),5),n,round(float(bb),4)))
    tot=int((m&v).sum()); stats[d]=dict(aoi_px=int(m.sum()),covered=tot,retrieved=retrieved,retrieved_frac=round(retrieved/max(tot,1),3))
    print(d,stats[d],{k:int(P[k].sum()) for k in fam})
    dates.append(d); 
    for k in fam: pres.setdefault(k,[]).append(P[k])
out={}
for k in fam:
    S=np.stack(pres[k]); cnt=S.sum(0); p=S.reshape(len(dates),-1).mean(1)
    probs=np.zeros(len(dates)+1); probs[0]=1
    for pk in p: probs=np.r_[probs[:1]*(1-pk),probs[1:]*(1-pk)+probs[:-1]*pk]
    out[k]=dict(per_date=[round(float(x),4) for x in p],obs_ge2=int((cnt>=2).sum()),exp_ge2=round(float(probs[2:].sum()*cnt.size),1),obs_ge3=int((cnt>=3).sum()),exp_ge3=round(float(probs[3:].sum()*cnt.size),2))
    print(k,out[k])
np.savez_compressed('data/emit_min_persist.npz',**{k:np.stack(pres[k]).sum(0).astype('uint8') for k in fam},dates=np.array(dates),grid=np.array([AOI[0],AOI[3],RES]))
json.dump(dict(stats=stats,persist=out,hits=hits),open('emit_min.json','w'))
c=collections.Counter(h[4] for h in hits); print(c.most_common(20))
