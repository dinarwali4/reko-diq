import numpy as np, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pyproj import Transformer
z=np.load('data/s2_change.npz'); B=z['B']; E=z['rgb_early']; L=z['rgb_late']; C=json.load(open('s2_change.json'))[:6]
fw=Transformer.from_crs(4326,32641,always_xy=True)
fig,axs=plt.subplots(2,6,figsize=(18,6.2),dpi=110)
g=lambda a: np.nan_to_num(np.clip(np.transpose(a,(1,2,0))/0.45,0,1)**(1/1.3))
for j,c in enumerate(C):
    x,y=fw.transform(c['centroid'][1],c['centroid'][0]); r=int((B[3]-y)/20); cc=int((x-B[0])/20); h=60
    sl=(slice(max(r-h,0),r+h),slice(max(cc-h,0),cc+h))
    axs[0,j].imshow(g(E[:,sl[0],sl[1]])); axs[1,j].imshow(g(L[:,sl[0],sl[1]]))
    axs[0,j].set_title(f"{c['centroid'][0]:.4f}N {c['centroid'][1]:.4f}E\n{c['ha']} ha",fontsize=8)
    for a in axs[:,j]: a.set_xticks([]); a.set_yticks([])
axs[0,0].set_ylabel('Apr–Aug 2025'); axs[1,0].set_ylabel('Jun–Sep 2026')
plt.tight_layout(); plt.savefig('fig/change_chips.png')
