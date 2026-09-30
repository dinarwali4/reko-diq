import numpy as np, rasterio, pystac_client, planetary_computer as pc
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from pyproj import Transformer
BB=[62.30,28.88,62.55,29.06]
tr=Transformer.from_crs(4326,32641,always_xy=True)
xs,ys=tr.transform([BB[0],BB[2],BB[0],BB[2]],[BB[1],BB[1],BB[3],BB[3]])
x0,y0,x1,y1=min(xs)//30*30,min(ys)//30*30,max(xs)//30*30+30,max(ys)//30*30+30
W,H=int((x1-x0)/30),int((y1-y0)/30); T=from_origin(x0,y1,30,30)
c=pystac_client.Client.open('https://planetarycomputer.microsoft.com/api/stac/v1',modifier=pc.sign_inplace)
out={}
for yr in [1990,1995,2000,2003,2006,2009,2011,2014,2017,2020,2022,2024,2025,2026]:
    its=[i for i in c.search(collections=['landsat-c2-l2'],bbox=BB,datetime=f'{yr}-01-01/{yr}-12-31',query={'eo:cloud_cover':{'lt':10}}).items()
         if i.properties.get('platform') in ('landsat-5','landsat-7','landsat-8','landsat-9')]
    if not its: print(yr,'none'); continue
    # prefer L5/L8/L9 over L7 (SLC-off gaps after 2003)
    its=sorted(its,key=lambda i:(i.properties['platform']=='landsat-7' and yr>2003, i.properties['eo:cloud_cover']))
    it=its[0]; bands={}
    for nm in ('red','green','blue','nir08','swir16','swir22'):
        if nm not in it.assets: continue
        with rasterio.open(it.assets[nm].href) as s, WarpedVRT(s,crs='EPSG:32641',transform=T,width=W,height=H,resampling=Resampling.bilinear) as v:
            a=v.read(1).astype('float32')
        bands[nm]=np.where(a>0,a*0.0000275-0.2,np.nan)
    out[yr]=np.stack([bands[k] for k in ('red','green','blue','nir08','swir16','swir22')])
    print(yr,it.id,it.properties['eo:cloud_cover'],flush=True)
np.savez_compressed('data/landsat.npz',years=np.array(list(out)),stack=np.stack(list(out.values())),T=np.array(list(T)[:6]))
print('done')
