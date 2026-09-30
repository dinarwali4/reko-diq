import numpy as np, json
from skimage import measure
from scipy import ndimage as ndi
from shapely.geometry import Polygon
F=np.load('data/facies.npz'); glon,glat=F['glon'],F['glat']; Z=json.load(open('facies_zones.json'))
sty={'phy':('ff2e9cff','552e9cff','Phyllic (sericite)'),'arg':('ffd93ad9','66d93ad9','Argillic (kaolinite/alunite)'),'cap':('ff1ad9f2','551ad9f2','Leached cap (jarosite + Fe)')}
k=['<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Reko Diq — RegenX alteration facies (Sep 2026)</name>']
for s,(l,p,_) in sty.items(): k.append(f'<Style id="{s}"><LineStyle><color>{l}</color><width>2.5</width></LineStyle><PolyStyle><color>{p}</color></PolyStyle></Style>')
k.append('<Style id="lab"><IconStyle><scale>0.6</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon></IconStyle><LabelStyle><color>ffffffff</color><scale>1.0</scale></LabelStyle></Style>')
k.append('<Style id="giv"><IconStyle><color>ffffe500</color><scale>1.2</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/wht-stars.png</href></Icon></IconStyle><LabelStyle><color>ffffe500</color></LabelStyle></Style>')
k.append('<Placemark><name>Given coordinate 28.972036, 62.432438</name><styleUrl>#giv</styleUrl><Point><coordinates>62.432438,28.972036,0</coordinates></Point></Placemark>')
for s,(l,p,name) in sty.items():
    k.append(f'<Folder><name>{name}</name>')
    m=ndi.binary_closing(ndi.binary_dilation(F[s],iterations=1),iterations=1)
    for c in measure.find_contours(np.pad(m.astype(float),1),0.5):
        poly=Polygon([(glon[0]+(cc-1)*0.0005,glat[0]-(r-1)*0.0005) for r,cc in c]).simplify(0.0002)
        if not poly.is_valid or poly.area<1.5e-7: continue
        co=' '.join(f'{x:.6f},{y:.6f},0' for x,y in poly.exterior.coords)
        k.append(f'<Placemark><styleUrl>#{s}</styleUrl><Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>{co}</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>')
    k.append('</Folder>')
k.append('<Folder><name>Zone labels (≥10 ha)</name>')
for z in Z:
    if z['area_ha']<10: continue
    k.append(f'<Placemark><name>{z["id"]}</name><description>{z["facies"]}; {z["area_ha"]} ha; centroid {z["centroid"][0]}N {z["centroid"][1]}E; {z["dist_given_km"]} km from given point.</description><styleUrl>#lab</styleUrl><Point><coordinates>{z["centroid"][1]},{z["centroid"][0]},0</coordinates></Point></Placemark>')
k.append('</Folder><Folder><name>Construction 2025–26 (Sentinel-2 change)</name>')
for n,la,lo in (('Airstrip',28.918,62.542),('Camp expansion',29.0421,62.5234),('New ponds / pads',29.0502,62.4897),('New corridor',28.8693,62.3113)):
    k.append(f'<Placemark><name>{n}</name><styleUrl>#lab</styleUrl><Point><coordinates>{lo},{la},0</coordinates></Point></Placemark>')
k.append('</Folder></Document></kml>')
open('RekoDiq_RegenX_alteration.kml','w').write('\n'.join(k)); print('ok')
