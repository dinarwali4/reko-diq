import numpy as np, netCDF4, glob, os
AOI=(62.30,28.88,62.55,29.06); RES=0.0005
NX=int(round((AOI[2]-AOI[0])/RES)); NY=int(round((AOI[3]-AOI[1])/RES))
glon=AOI[0]+(np.arange(NX)+.5)*RES; glat=AOI[3]-(np.arange(NY)+.5)*RES
for f in sorted(glob.glob('data/emit/EMIT_L2A_RFL_001_*.nc')):
    date=f.split('_')[4][:8]; out=f'data/emit_aoi_{date}.npz'
    if os.path.exists(out): continue
    try: ds=netCDF4.Dataset(f); ds.set_auto_mask(False)
    except OSError: print(date,'incomplete'); continue
    wl=ds['sensor_band_parameters/wavelengths'][:]; gx=ds['location/glt_x'][:]; gy=ds['location/glt_y'][:]; gt=ds.getncattr('geotransform')
    oj=((glon-gt[0])/gt[1]).astype(int); oi=((glat-gt[3])/gt[5]).astype(int); OI,OJ=np.meshgrid(oi,oj,indexing='ij')
    ins=(OI>=0)&(OJ>=0)&(OI<gx.shape[0])&(OJ<gx.shape[1]); X=np.zeros((NY,NX),int); Y=np.zeros((NY,NX),int)
    X[ins]=gx[OI[ins],OJ[ins]]; Y[ins]=gy[OI[ins],OJ[ins]]; v=ins&(X>0)&(Y>0)
    print(date,'covered',round(float(v.mean()),3),flush=True)
    if v.mean()<0.2: continue
    rows=Y[v]-1; cols=X[v]-1; r0,r1=rows.min(),rows.max()+1
    R=ds['reflectance'][r0:r1,:,:]; cube=np.full((NY,NX,R.shape[2]),np.nan,'float32'); cube[v]=R[rows-r0,cols]; cube[cube<-9000]=np.nan
    m=netCDF4.Dataset(f.replace('_RFL_','_MASK_')); m.set_auto_mask(False); agg=m['mask'][:,:,7]; mk=np.ones((NY,NX)); mk[v]=agg[rows,cols]
    np.savez_compressed(out,cube=cube.astype('float16'),wl=wl,mask=mk.astype('uint8'),glon=glon,glat=glat); print(date,'saved',flush=True)
