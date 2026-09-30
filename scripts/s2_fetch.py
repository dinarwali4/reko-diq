import numpy as np, rasterio, pystac_client, concurrent.futures as cf
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from pyproj import Transformer
BBOX=[62.25,28.85,62.60,29.10]
c=pystac_client.Client.open('https://earth-search.aws.element84.com/v1')
items=[i for i in c.search(collections=['sentinel-2-l2a'],bbox=BBOX,datetime='2024-10-01/2026-09-25',query={'eo:cloud_cover':{'lt':2}}).items()
       if True]
tiles={}
for i in items: tiles.setdefault(i.id.split('_')[1],[]).append(i)
print({k:len(v) for k,v in tiles.items()})
tile=max(tiles,key=lambda k:len(tiles[k])); its=sorted(tiles[tile],key=lambda i:i.properties['eo:cloud_cover'])[:14]
print('tile',tile,[i.id for i in its])
epsg=32641; tr=Transformer.from_crs(4326,epsg,always_xy=True)
xs,ys=tr.transform([BBOX[0],BBOX[2],BBOX[0],BBOX[2]],[BBOX[1],BBOX[1],BBOX[3],BBOX[3]])
B=(min(xs)//20*20,min(ys)//20*20,max(xs)//20*20+20,max(ys)//20*20+20); W=int((B[2]-B[0])/20); H=int((B[3]-B[1])/20)
bands={'blue':'B02','green':'B03','red':'B04','rededge1':'B05','rededge2':'B06','rededge3':'B07','nir':'B08','nir08':'B8A','swir16':'B11','swir22':'B12','scl':'SCL'}
def rd(href,res):
    with rasterio.open(href) as s:
        return s.read(1,window=from_bounds(*B,s.transform),out_shape=(H,W),resampling=res,boundless=True,fill_value=0)
def job(it): return it.id,np.stack([rd(it.assets[k].href,Resampling.nearest if k=='scl' else Resampling.average) for k in bands])
out={}
with cf.ThreadPoolExecutor(6) as ex:
    for iid,a in ex.map(job,its): out[iid]=a; print('got',iid,flush=True)
np.savez_compressed('data/s2_raw.npz',ids=np.array(list(out)),stack=np.stack(list(out.values())),bounds=np.array(B),epsg=epsg,bandnames=np.array(list(bands.values())),
   boa=np.array([i.properties.get('earthsearch:boa_offset_applied') for i in its]))
print('done',W,H,epsg)
