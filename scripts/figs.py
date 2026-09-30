import numpy as np, json, rasterio, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from matplotlib.patches import Patch
z=np.load('data/emit_aoi_20231001.npz'); C=z['cube']; wl=z['wl']; glon,glat=z['glon'],z['glat']; bi=lambda w:int(np.argmin(abs(wl-w)))
rgb=np.nan_to_num(np.clip(np.stack([C[...,bi(650)],C[...,bi(560)],C[...,bi(470)]],-1).astype(float)/0.45,0,1)**(1/1.3))
G=np.load('data/zone_grid.npz'); D=G['D'].astype(float); gy,gx=np.gradient(D,55.3,48.6); sl=np.arctan(np.hypot(gx,gy)); asp=np.arctan2(-gx,gy)
hs=np.clip(np.sin(np.radians(45))*np.cos(sl)+np.cos(np.radians(45))*np.sin(sl)*np.cos(np.radians(315)-asp),0,1)
F=np.load('data/facies.npz'); ext=[glon[0],glon[-1],glat[-1],glat[0]]; asp_=1/np.cos(np.radians(28.97))
fig,ax=plt.subplots(figsize=(13,9.6),dpi=130); ax.imshow(rgb*(0.6+0.4*hs[...,None]),extent=ext,aspect=asp_)
cols={'phy':(1.0,0.62,0.18,0.75),'arg':(0.85,0.2,0.85,0.85),'cap':(0.95,0.85,0.1,0.8)}
for k in ('phy','cap','arg'):
    ov=np.zeros(F[k].shape+(4,)); ov[F[k]]=cols[k]; ax.imshow(ov,extent=ext,aspect=asp_,interpolation='nearest')
ax.plot(62.432438,28.972036,'*',ms=18,color='#00e5ff',mec='k'); ax.text(62.435,28.966,'Given coordinate\n28.972036, 62.432438',fontsize=8,color='w',path_effects=[pe.withStroke(linewidth=2,foreground='k')])
feat=[('Airstrip (new, 2025–26)',28.918,62.542,'s'),('Camp expansion',29.0421,62.5234,'^'),('New ponds/pads',29.0502,62.4897,'D'),('New corridor',28.8693,62.3113,'v')]
for n,la,lo,m in feat:
    ax.plot(lo,la,m,ms=10,color='w',mec='k'); t=ax.text(lo+(0.004 if lo<62.53 else -0.004),la,n,ha=('left' if lo<62.53 else 'right'),fontsize=7.5,color='w',path_effects=[pe.withStroke(linewidth=2,foreground='k')]); t.set_clip_on(True)
for n,(la,lo) in {'Western complex':(28.985,62.385),'North-west cap':(29.033,62.39),'Far-west phyllic':(29.045,62.32)}.items():
    ax.text(lo,la,n,fontsize=9,color='w',weight='bold',ha='center',path_effects=[pe.withStroke(linewidth=2.5,foreground='k')])
ax.add_patch(plt.Rectangle((62.30,28.88),0.10,0.04,fill=False,ec='w',ls='--',lw=1)); ax.text(62.302,28.883,'Alluvial control: 0 zone cells',fontsize=7.5,color='w',path_effects=[pe.withStroke(linewidth=2,foreground='k')])
ax.legend(handles=[Patch(color=cols['phy'],label='Phyllic (sericite, 2200 nm)'),Patch(color=cols['arg'],label='Argillic (kaolinite/alunite, 2165 nm)'),Patch(color=cols['cap'],label='Leached cap (jarosite 2265 nm + Fe³⁺)')],loc='lower right',fontsize=8,framealpha=.9)
ax.set_xlim(ext[0],ext[1]); ax.set_ylim(ext[2],ext[3]); ax.set_xlabel('Longitude °E'); ax.set_ylabel('Latitude °N')
ax.set_title('Reko Diq: porphyry alteration facies persistent on both Oct 2023 and Oct 2024 EMIT dates',fontsize=10)
plt.tight_layout(); plt.savefig('fig/facies_map.png',bbox_inches='tight')
# spectra
fig,axs=plt.subplots(1,2,figsize=(13,4.6),dpi=130)
LON,LAT=np.meshgrid(glon,glat); ctrl=(LON>=62.30)&(LON<=62.40)&(LAT>=28.88)&(LAT<=28.92)
for ax,d in zip(axs,['20231001','20250202']):
    cc=np.load(f'data/emit_aoi_{d}.npz'); Cc=cc['cube']; w=cc['wl']; sel=(w>2100)&(w<2400); ww=w[sel]
    for k,(n,m,col) in enumerate((('Argillic zones',F['arg'],'#b03ab0'),('Phyllic zones',F['phy'],'#e08a2e'),('Leached-cap zones',F['cap'],'#c9a800'),('Alluvial control',ctrl,'#6c7a74'))):
        s=np.nanmean(Cc[m][:,sel].astype('float32'),0); ct=s[0]+(s[-1]-s[0])*(ww-ww[0])/(ww[-1]-ww[0]); ax.plot(ww,s/ct-0.03*k,color=col,lw=1.8,label=f'{n} (n={int(m.sum())})')
    for x,t in ((2165,'2165 kaol/alunite'),(2200,'2200 sericite'),(2265,'2265 jarosite'),(2330,'2330 carbonate/chlorite')): ax.axvline(x,color='#999',ls='--',lw=.7); ax.text(x+1,ax.get_ylim()[1],t,fontsize=7,rotation=90,va='top')
    ax.set_title(f'EMIT {d[:4]}-{d[4:6]}-{d[6:]}'+(' (independent test date)' if d.startswith('2025') else ''),fontsize=10); ax.set_xlabel('Wavelength (nm)'); ax.spines[['top','right']].set_visible(False); ax.grid(alpha=.2)
axs[0].set_ylabel('Continuum-removed (offset 0.03)'); axs[1].legend(frameon=False,fontsize=7.5,loc='lower left')
plt.tight_layout(); plt.savefig('fig/spectra.png'); print('ok')
