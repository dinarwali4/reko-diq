import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
z=np.load('data/emit_feat.npz'); glon,glat=z['glon'],z['glat']
def rz(a): med=np.nanmedian(a); mad=1.4826*np.nanmedian(abs(a-med)); return (a-med)/mad
ext=[glon[0],glon[-1],glat[-1],glat[0]]
fig,axs=plt.subplots(2,3,figsize=(16,9.5),dpi=110)
C=np.load('data/emit_aoi_20231001.npz'); cube=C['cube']; wl=C['wl']; bi=lambda w:int(np.argmin(abs(wl-w)))
rgb=np.nan_to_num(np.clip(np.stack([cube[...,bi(650)],cube[...,bi(560)],cube[...,bi(470)]],-1).astype(float)/0.45,0,1)**(1/1.3))
axs[0,0].imshow(rgb,extent=ext); axs[0,0].set_title('EMIT true colour, 1 Oct 2023')
for ax,(k,t) in zip(axs.ravel()[1:],[('AlOH_2200','Al-OH 2200 (sericite)'),('Kaol_2165','2165 (kaolinite/alunite)'),('Jar_2265','Jarosite 2265'),('Fe3_900','Fe³⁺ 900'),('MgOH_2330','Mg-OH/carbonate 2330')]):
    Zm=np.nanmean([rz(z[f'{d}_{k}']) for d in ('20231001','20241001')],0)
    im=ax.imshow(Zm,extent=ext,cmap='magma',vmin=0,vmax=4); ax.set_title(t+' — mean z, Oct 2023 & Oct 2024'); plt.colorbar(im,ax=ax,fraction=.035)
for ax in axs.ravel(): ax.plot(62.432438,28.972036,'c*',ms=12,mec='k')
plt.tight_layout(); plt.savefig('fig/quicklook.png')
