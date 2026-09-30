import json
from PIL import Image, ImageDraw, ImageFont
F=lambda s,b=True: ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial Bold.ttf' if b else '/System/Library/Fonts/Supplemental/Arial.ttf',s)
Z={z['id']:z for z in json.load(open('facies_zones.json'))}
jobs=[('RekoDiq_A1_facies_lease_2023-09','Reko Diq: persistent porphyry alteration facies','Airbus imagery Sep 2023 · EMIT hyperspectral Oct 2023 + Oct 2024, confirmed Feb 2025',['RZ1','RZ2','RZ3','RZ4','RZ5']),
 ('RekoDiq_A2_facies_west_2023-09','Western porphyry cluster: phyllic, argillic and leached-cap zones','Airbus Sep 2023 · zones rendered by Google Earth on terrain',['RZ2','RZ3','RZ4','RZ7','RZ8']),
 ('RekoDiq_H1_west_2003-12','Dec 2003 · before the main drilling campaigns','Same frame and zones in all historical views',[]),
 ('RekoDiq_H2_west_2010-08','Aug 2010 · Tethyan Copper feasibility-era drilling','Same frame and zones in all historical views',[]),
 ('RekoDiq_H3_west_2026-09','Sep 2026 · latest Google Earth mosaic','Same frame and zones in all historical views',[]),
 ('RekoDiq_H1_center_2004-04','Given coordinate, Apr 2004','28.972036°N 62.432438°E',[]),
 ('RekoDiq_H2_center_2010-11','Given coordinate, Nov 2010','28.972036°N 62.432438°E',[])]
for f,t,s,ids in jobs:
    im=Image.open(f'ge/img/{f}.jpg').convert('RGB'); w,h=im.size; d=ImageDraw.Draw(im)
    d.rectangle([0,0,w,122],fill=(24,20,16)); d.text((28,14),t,font=F(40),fill=(255,255,255)); d.text((30,66),s,font=F(24,False),fill=(214,206,196))
    d.text((w-560,24),'RegenX.eco Geospatial Computational Model',font=F(22),fill=(143,193,83)); d.text((w-560,58),'Prepared by Dinar Wali · Sep 2026',font=F(22,False),fill=(214,206,196))
    if 'center' not in f:
        d.rectangle([w-470,122,w,238],fill=(24,20,16))
        for i,(c,l) in enumerate((((255,156,46),'Phyllic (sericite)'),((217,58,217),'Argillic (kaolinite/alunite)'),((242,217,26),'Leached cap (jarosite + Fe)'))):
            d.rectangle([w-450,140+i*32,w-420,158+i*32],outline=c,width=4); d.text((w-408,136+i*32),l,font=F(20,False),fill=(235,235,235))
    if ids:
        bx,by=24,h-60-44*len(ids)-50; d.rectangle([bx-10,by-10,bx+700,h-40],fill=(24,20,16))
        d.text((bx,by),'Zone centroids (WGS84, ±60 m) · facies · area',font=F(22),fill=(255,255,255))
        for i,k in enumerate(ids):
            z=Z[k]; d.text((bx,by+40+i*44),f"{k:5s} {z['centroid'][0]:.5f}°N {z['centroid'][1]:.5f}°E · {z['facies'].split(' (')[0]} · {z['area_ha']} ha",font=F(23,False),fill=(255,190,110))
    im.save(f'fig/{f.replace("RekoDiq_","annot_")}.jpg',quality=90)
# history panel
tiles=[Image.open(f'fig/annot_{n}.jpg').crop((0,122,2398,1522)).resize((999,583)) for n in ('H1_west_2003-12','H2_west_2010-08','H3_west_2026-09','A2_facies_west_2023-09')]
labs=['Dec 2003','Aug 2010','Sep 2026 mosaic','Sep 2023 + zones']
W,H=999,583; out=Image.new('RGB',(W*2+12,H*2+12+140),(24,20,16)); dd=ImageDraw.Draw(out)
dd.text((20,18),'Reko Diq western porphyries, 2003 → 2026 (Google Earth historical imagery)',font=F(38),fill=(255,255,255))
dd.text((20,72),'The alteration zones are drawn on every epoch to fix position. RegenX.eco · Prepared by Dinar Wali',font=F(24,False),fill=(214,206,196))
for i,(t,l) in enumerate(zip(tiles,labs)):
    x,y=(i%2)*(W+12),140+(i//2)*(H+12); out.paste(t,(x,y)); ImageDraw.Draw(out).rectangle([x,y,x+330,y+50],fill=(24,20,16)); ImageDraw.Draw(out).text((x+12,y+8),l,font=F(30),fill=(255,255,255))
out.save('fig/history_panel.jpg',quality=88); print('ok')
