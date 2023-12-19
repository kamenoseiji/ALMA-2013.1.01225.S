prefix =  'uid___A002_Xa76868_X2423'
REFANT = 'DA55'
BPCAL1 = 'J0423-0120'
BPCAL2 = 'J0238+166'
PHCAL = 'J0239-0234'
TARGET = 'NGC_1052'
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0423-0120          04:23:15.800730 -01.20.33.06550 J2000   0        2546081
#  1    none J0238+166           02:38:38.930107 +16.36.59.27462 J2000   1        1197949
#  2    none NGC_1052            02:41:04.798500 -08.15.20.75180 J2000   2         637510
#  3    none J0239-0234          02:39:45.472270 -02.34.40.91440 J2000   3         525807
#  4    none NGC_1052            02:41:04.798510 -08.15.20.75170 J2000   4        6008652
os.system('rm -rf ' + prefix + '_flagonline.txt')    
importasdm(prefix)
browsetable(tablename = prefix+'.ms')
listobs(vis=prefix+'.ms', scan='', spw='', verbose=True, listfile=prefix+'.listobs')   # SPW=[17,19,21,23]
plotants( vis=prefix + '.ms', figfile=prefix + '_plotants.png') #REFANT = DA55
#
#
flagdata(vis=prefix + '.ms', mode='manual', spw='0~24', autocorr=True, flagbackup = False)
flagdata(vis=prefix + '.ms', mode='manual', intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*', flagbackup = False)
flagdata(vis=prefix + '.ms', mode='manual', timerange='2015/08/05/08:28:52~08:29:00', flagbackup=False)
flagdata(vis=prefix + '.ms', mode='manual', timerange='2015/08/05/08:39:52~08:39:54', flagbackup=False)
flagdata(vis=prefix + '.ms', mode='manual', timerange='2015/08/05/08:56:58~08:57:07', flagbackup=False)
flagdata(vis=prefix + '.ms', mode='manual', timerange='2015/08/05/09:01:52.5~09:01:56', flagbackup=False)
#flagdata(vis=prefix + '.ms', mode='summary', name='after')
flagmanager(vis=prefix+'.ms', mode='save', versionname='Apriori')
#
#-------- Tsys 
from recipes.almahelpers import tsysspwmap
os.system('rm -rf cal.' + prefix + '.tsys')
gencal(vis=prefix+'.ms', caltype='tsys', caltable='cal.'+prefix+'.tsys')
plotcal(caltable='cal.'+prefix+'.tsys', xaxis='freq', yaxis='tsys', spw='9:3~124,11:3~124,13:3~124,15:3~124', iteration='antenna', subplot=221)
#-------- WVR 
os.system('rm -rf cal.' + prefix + '.WVR')
wvrgcal(segsource=True, caltable='cal.'+prefix+'.WVR', vis=prefix+'.ms', wvrspw=[0], tie=[BPCAL1 + ',' + BPCAL2+ ',' + PHCAL+ ',' + TARGET], toffset=0, wvrflag=[], statsource=TARGET)
plotcal('cal.'+prefix+'.WVR', xaxis='time',yaxis='phase',iteration='antenna', subplot=221)
#
tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable='cal.'+prefix+'.tsys', tsysChanTol=1)
for fields in [BPCAL1, BPCAL2, PHCAL, TARGET]:
    applycal(vis=prefix+'.ms', field=fields, flagbackup=False, spw='17,19,21,23', interp=['linear', 'nearest'], gaintable=['cal.'+prefix+'.tsys', 'cal.'+prefix+'.WVR'], gainfield=fields, spwmap=[tsysmap, []], calwt=True)
#-------- Split
os.system('rm -rf '+prefix+'_line.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_line.ms', datacolumn='corrected', spw='17,19,21,23')
#
#-------- Phase Cal for bandpass
os.system('rm -rf P0*')
gaincal(vis=prefix+'_line.ms', caltable='P0', spw='*:4~123', field= BPCAL2, scan='', selectdata=True, solint='int', refant=REFANT, calmode='p')
plotcal(caltable = 'P0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-180,180], iteration = 'antenna', figfile='P0.phase.png', subplot = 221)
os.system('rm -rf B0*')
bandpass(vis = prefix + '_line.ms', caltable = 'B0', gaintable='P0', spw='*:3~124', field= BPCAL2, scan='', minblperant=5, minsnr=5, solint='inf', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True, spwmap=[0,1,2,3])
plotbandpass(caltable = 'B0', xaxis='freq', yaxis='amp', plotrange = [0,0,0,1.2], figfile='B0.png')
plotbandpass(caltable = 'B0', xaxis='freq', yaxis='phase',plotrange = [0,0,-180.0,180.0], figfile='B0.png')
#-------- Phase Cal for all
os.system('rm -rf P1*')
gaincal(vis=prefix+'_line.ms', caltable='P1', spw='*:3~124', gaintable = ['B0','P0'], field=BPCAL1 + ', ' +  BPCAL2 + ', ' + PHCAL + ', ' + TARGET, selectdata=True, solint='int', refant=REFANT, gaintype='T', combine='spw', calmode='p', minsnr=3)
plotcal(caltable = 'P1', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-180,180], iteration = 'antenna', figfile='cal_phase.png', subplot = 221)
#-------- Flux cal
au.getALMAFlux(BPCAL2, 230.07939217, date='2015-08-15', showplot=True)
#au.getALMAFlux(BPCAL2, 230.07939217, date='2015-08-15', showplot=True)
#Closest Band 3 measurement: 4.080 +- 0.200 (age=+2 days) 103.5 GHz
#Closest Band 3 measurement: 4.200 +- 0.140 (age=+2 days) 91.5 GHz
#Closest Band 7 measurement: 2.350 +- 0.080 (age=-1 days) 337.5 GHz
#getALMAFluxCSV(): Fitting for spectral index with 1 measurement pair of age -12 days from 2015-08-15, with age separation of 0 days
#  2015-08-27: freqs=[91.46, 103.49, 343.48], fluxes=[3.53, 3.42, 1.88]
#Median Monte-Carlo result for 230.079392 = 2.284187 +- 0.156257 (scaled MAD = 0.156264)
#Result using spectral index of -0.482145 for 230.079 GHz from 4.140 Jy at 97.475 GHz = 2.736325 +- 0.156257 Jy
os.system('rm -rf G0*')
setjy( vis = prefix + '_line.ms', field=BPCAL2, spw='0,1,2,3', standard='manual', fluxdensity=[2.736325,0,0,0], spix = [-0.482,0], reffreq = '230.07939217GHz', usescratch=False)
gaincal(vis = prefix + '_line.ms', caltable = 'G0', spw ='*:3~124', field = '', minsnr=5.0, solint='inf', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['B0','P0','P1'], spwmap=[[0,1,2,3],[0,1,2,3],[0,0,0,0]], calmode = 'a')
plotcal(caltable = 'G0', xaxis = 'time', yaxis = 'amp', plotsymbol='o', plotrange = [], iteration = 'antenna', figfile='G0.png', subplot = 221)
fluxscale(vis= prefix + '_line.ms', caltable='G0', fluxtable='G0.flux', reference=BPCAL2, transfer=BPCAL1 + ', ' +  PHCAL + ', ' + TARGET, refspwmap=[0,1,2,3])
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0423-0120 in SpW=0 (freq=2.29428e+11 Hz) is: 1.13398 +/- 0.00226296 (SNR = 501.106, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0423-0120 in SpW=1 (freq=2.30794e+11 Hz) is: 1.13141 +/- 0.00214056 (SNR = 528.556, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0423-0120 in SpW=2 (freq=2.1606e+11 Hz) is: 1.19524 +/- 0.00228794 (SNR = 522.41, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0423-0120 in SpW=3 (freq=2.14064e+11 Hz) is: 1.19863 +/- 0.00251143 (SNR = 477.27, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=0 (freq=2.29428e+11 Hz) is: 0.477146 +/- 0.000823216 (SNR = 579.613, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=1 (freq=2.30794e+11 Hz) is: 0.477251 +/- 0.000846724 (SNR = 563.644, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=2 (freq=2.1606e+11 Hz) is: 0.51101 +/- 0.000819667 (SNR = 623.436, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=3 (freq=2.14064e+11 Hz) is: 0.514843 +/- 0.000807545 (SNR = 637.54, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0239-0234 in SpW=0 (freq=2.29428e+11 Hz) is: 0.182447 +/- 0.000941478 (SNR = 193.788, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0239-0234 in SpW=1 (freq=2.30794e+11 Hz) is: 0.180728 +/- 0.000864768 (SNR = 208.991, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0239-0234 in SpW=2 (freq=2.1606e+11 Hz) is: 0.191157 +/- 0.000899997 (SNR = 212.397, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for J0239-0234 in SpW=3 (freq=2.14064e+11 Hz) is: 0.192492 +/- 0.000920575 (SNR = 209.1, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=0 (freq=2.29428e+11 Hz) is: 0.475096 +/- 0.000768366 (SNR = 618.32, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=1 (freq=2.30794e+11 Hz) is: 0.475103 +/- 0.000743445 (SNR = 639.056, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=2 (freq=2.1606e+11 Hz) is: 0.509164 +/- 0.000759129 (SNR = 670.722, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Flux density for NGC_1052 in SpW=3 (freq=2.14064e+11 Hz) is: 0.513157 +/- 0.000754785 (SNR = 679.871, N = 74)
#2019-09-23 00:05:24 INFO fluxscale	 Fitted spectrum for J0423-0120 with fitorder=1: Flux density = 1.16453 +/- 0.00154912 (freq=222.458 GHz) spidx: a_1 (spectral index) =-0.811383 +/- 0.0392169 covariance matrix for the fit:  covar(0,0)=9.71106e-07 covar(0,1)=-2.45074e-06 covar(1,0)=-2.45074e-06 covar(1,1)=0.00447483
#2019-09-23 00:05:24 INFO fluxscale	 Fitted spectrum for NGC_1052 with fitorder=1: Flux density = 0.49472 +/- 0.000833282 (freq=222.458 GHz) spidx: a_1 (spectral index) =-1.06171 +/- 0.0494341 covariance matrix for the fit:  covar(0,0)=6.97169e-07 covar(0,1)=4.69745e-06 covar(1,0)=4.69745e-06 covar(1,1)=0.00318388
#2019-09-23 00:05:24 INFO fluxscale	 Fitted spectrum for J0239-0234 with fitorder=1: Flux density = 0.186618 +/- 0.000212655 (freq=222.458 GHz) spidx: a_1 (spectral index) =-0.815699 +/- 0.0334007 covariance matrix for the fit:  covar(0,0)=5.88952e-06 covar(0,1)=1.59327e-05 covar(1,0)=1.59327e-05 covar(1,1)=0.0268272
#2019-09-23 00:05:24 INFO fluxscale	 Fitted spectrum for NGC_1052 with fitorder=1: Flux density = 0.492828 +/- 0.000811201 (freq=222.458 GHz) spidx: a_1 (spectral index) =-1.07349 +/- 0.0482287 covariance matrix for the fit:  covar(0,0)=5.90182e-07 covar(0,1)=2.79889e-06 covar(1,0)=2.79889e-06 covar(1,1)=0.00268635
#2019-09-23 00:05:24 INFO fluxscale	Storing result in G0.flux
#2019-09-23 00:05:24 INFO fluxscale	Writing solutions to table: G0.flux
#
#
applycal(vis= prefix + '_line.ms', flagbackup=False, field='', interp=['nearest','nearest', 'linear','nearest'], gainfield = '', gaintable=['B0','P0','P1', '.G0.flux'], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]])
#-------- Split into target source
os.system('rm -rf ' + BPCAL1 + '.ms'); split(vis=prefix + '_line.ms', outputvis=BPCAL1+'.ms', field=BPCAL1, datacolumn='corrected')
os.system('rm -rf ' + BPCAL2 + '.ms'); split(vis=prefix + '_line.ms', outputvis=BPCAL2+'.ms', field=BPCAL2, datacolumn='corrected')
os.system('rm -rf ' + PHCAL + '.ms'); split(vis=prefix + '_line.ms', outputvis=PHCAL+'.ms', field=PHCAL, datacolumn='corrected')
os.system('rm -rf ' + TARGET + '.ms'); split(vis=prefix + '_line.ms', outputvis=TARGET+'.ms', field=TARGET, datacolumn='corrected')
#-------- Imaging calibrators
##---- Imaging J0423-0120
os.system('rm -rf SC.Gap0')
gaincal(vis=BPCAL1+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[1.12920191712223, 0.0224139593914416, -0.0439291237833308, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=BPCAL1+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + BPCAL1 + '_SPW*')
tclean( vis=BPCAL1+'.ms', imagename=BPCAL1 + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL1+'.ms', imagename=BPCAL1 + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL1+'.ms', imagename=BPCAL1 + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL1+'.ms', imagename=BPCAL1 + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
##---- Imaging J0238+166
os.system('rm -rf SC.Gap0')
gaincal(vis=BPCAL2+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[2.69927, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=BPCAL2+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + BPCAL2 + '_SPW*')
tclean( vis=BPCAL2+'.ms', imagename=BPCAL2 + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL2+'.ms', imagename=BPCAL2 + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL2+'.ms', imagename=BPCAL2 + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=BPCAL2+'.ms', imagename=BPCAL2 + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
##---- Imaging J0239-0234
os.system('rm -rf SC.Gap0')
gaincal(vis=PHCAL+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[0.183757, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=PHCAL+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + PHCAL + '_SPW*')
tclean( vis=PHCAL+'.ms', imagename=PHCAL + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=PHCAL+'.ms', imagename=PHCAL + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=PHCAL+'.ms', imagename=PHCAL + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=PHCAL+'.ms', imagename=PHCAL + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
#
#-------- Imaging NGC1052
# Line-free CH range
# SPW0 : 0:2~6;32~50;110~124
# SPW1 : 1:50~92;118~124
# SPW2 : 2:4~10;90~124
# SPW3 : 3:4~124
os.system('rm -rf SC.Gap0')
gaincal(vis=TARGET+'.ms', caltable='SC.Gap0', spw='0:3~6;32~50;110~124, 1:50~92;118~124, 2:4~10;90~124, 3:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[0.492828, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=5)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
applycal(vis=TARGET+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + TARGET + '_SPW*')
tclean(vis=TARGET+'.ms', datacolumn='corrected', imagename=TARGET+'_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='230.5380GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # CO J=2-1
tclean(vis=TARGET+'.ms', datacolumn='corrected', imagename=TARGET+'_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='231.90093GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # H30alpha
tclean(vis=TARGET+'.ms', datacolumn='corrected', imagename=TARGET+'_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='217.10498GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # SiO J=5-4 v=0
tclean(vis=TARGET+'.ms', datacolumn='corrected', imagename=TARGET+'_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='215.220653GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # SO 5_5-4_4
#-------- Continuum imaging
os.system('rm -rf ' + TARGET + '_Cont*')
tclean(vis=TARGET+'.ms', datacolumn='corrected', imagename=TARGET+'_Cont', spw='0:3~6;32~50;110~124, 1:50~92;118~124, 2:4~10;90~124, 3:3~124', specmode='mfs', niter=5000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True)
imview(TARGET+'_Cont.image.pbcor')
#-------- Continuum subtraction
os.system('rm -rf ' + TARGET + 'ContSub.ms*')
uvcontsub(vis = TARGET+'.ms', spw='0', fitspw='0:3~6;32~50;110~124', solint ='int', fitorder = 1)
os.system('rm -rf ' + TARGET + '_ContSub*')
tclean(vis=TARGET+'.ms.contsub', datacolumn='corrected', imagename=TARGET+'_ContSub', spw='0', specmode='cube', start=34, nchan=45, width=1, outframe='LSRK', veltype='radio', restfreq='230.5380GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', uvtaper='1000klambda', threshold='0.000mJy', pbcor=True) # CO J=2-1
os.system('rm -rf ' + TARGET + '_ContSub.image.mom*')
immoments(TARGET+'_ContSub.image', moments=[0], outfile=TARGET + '_ContSub.image.mom0', chans='1~44', includepix=[3e-4, 100])
imview(raster={'file': TARGET + '_ContSub.image.mom0', 'range': [0.0, 0.2], 'colormap': 'RGB 2'}, contour={'file': TARGET+'_Cont.image.pbcor', 'levels': [1,2,4,8,16,32,64], 'unit': 0.0001}, zoom={'blc': [64,64], 'trc': [447,447]})
os.system('rm -rf ' + TARGET + '_ContSub.image.mom1')
immoments(TARGET+'_ContSub.image', moments=[1], outfile=TARGET + '_ContSub.image.mom1', chans='1~44', includepix=[7e-4, 100])
imview(raster={'file': TARGET + '_ContSub.image.mom1', 'range': [1200, 1900], 'colormap': 'Rainbow 3'}, contour={'file': TARGET+'_Cont.image.pbcor', 'levels': [1,2,4,8,16,32,64], 'unit': 0.0001}, zoom={'blc': [64,64], 'trc': [447,447]})
os.system('rm -rf ' + TARGET + '_ContSubforChanMap*')
tclean(vis=TARGET+'.ms.contsub', datacolumn='corrected', imagename=TARGET+'_ContSubforChanMap', spw='0', specmode='cube', start='1160km/s', nchan=16, width='45km/s', outframe='LSRK', veltype='radio', restfreq='230.5380GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', uvtaper='1000klambda', threshold='0.000mJy', pbcor=True) # CO J=2-1
imview(raster={'file': TARGET + '_ContSubforChanMap.image.pbcor', 'range': [0.000, 0.003], 'colormap': 'Rainbow 2'}, contour={'file': TARGET + '_ContSubforChanMap.image', 'levels': [-100, -10, 1,2,4,8,16,32,64,128,256,512,1024], 'unit': 0.00078}, zoom={'blc': [96,96], 'trc': [415,415]})
#-------- For Momoent-0 map
# image rms = 0.26 mJy/beam
os.system('rm -rf ' + TARGET + '_ContSubforChanMap.image.mom0')
immoments(TARGET+'_ContSubforChanMap.image', moments=[0], outfile=TARGET + '_ContSubforChanMap.image.mom0', chans='1~12', includepix=[0.00052, 100])
#imview(raster={'file': TARGET + '_ContSubforChanMap.image.mom0', 'range': [0.0, 0.25], 'colormap': 'Rainbow 2', 'colorwedge': True, 'title': 'NGC 1052 CO (J=2-1) + Continuum'}, zoom={'blc': [96,96], 'trc': [415,415]})
#imview(raster={'file': TARGET + '_ContSubforChanMap.image.mom0', 'range': [0.0, 0.25], 'colormap': 'Rainbow 2', 'colorwedge': True, 'title': 'NGC 1052 CO (J=2-1) + Continuum'}, contour={'file': TARGET+'_Cont.image.pbcor', 'levels': [1,2,4,8,16,32,64, 128, 256, 512, 1024, 2048, 4096], 'unit': 0.0001}, zoom={'blc': [96,96], 'trc': [415,415]})
#imview(raster={'file': TARGET + '_ContSubforChanMap.image.mom0', 'range': [0.0, 0.25], 'colormap': 'Rainbow 2', 'colorwedge': True, 'title': 'NGC 1052 CO (J=2-1) + Continuum'}, contour={'file': TARGET+'_Cont.image.pbcor', 'levels': [1,2,4,8,16,32,64, 128, 256, 512, 1024, 2048, 4096], 'unit': 0.0001}, zoom={'blc': [96,96], 'trc': [415,415]}, out={'file': 'NGC1052_CO21+cont.pdf', 'dpi': 300})
#-------- For Momoent-1 map
os.system('rm -rf ' + TARGET + '_ContSubforChanMap.image.mom1')
immoments(TARGET+'_ContSubforChanMap.image', moments=[1], outfile=TARGET + '_ContSubforChanMap.image.mom1', chans='1~12', includepix=[0.00078, 100])
imview(raster={'file': TARGET + '_ContSubforChanMap.image.mom1', 'range': [1200, 1800], 'colormap': 'Rainbow 3'}, contour={'file': TARGET+'_Cont.image.pbcor', 'levels': [1,2,4,8,16,32,64,128, 256, 512, 1024, 2048, 4096], 'unit': 0.0001}, zoom={'blc': [96,96], 'trc': [415,415]})
os.system('rm -rf ' + TARGET + '_ContSubforPV*')
tclean(vis=TARGET+'.ms.contsub', datacolumn='corrected', imagename=TARGET+'_ContSubforPV', spw='0', specmode='cube', start=34, nchan=45, width=1, outframe='LSRK', veltype='radio', restfreq='230.5380GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', uvtaper='1000klambda', threshold='0.000mJy', pbcor=True) # CO J=2-1
#-------- For PV map
imtrans(imagename="NGC_1052_ContSubforPV.image.pbcor",outfile="NGC_1052_ContSubforPVorder.image.pbcor",order="12-3")
imview(raster={'file': TARGET+'_ContSubforPVorder.image.pbcor', 'colormap': 'Hot metal 2'}, zoom={'blc': [96,96], 'trc': [415,415]})
imview(raster={'file': TARGET+'_ContSubforPVorder.image.pbcor', 'colormap': 'Hot metal 2'}, zoom={'blc': [96,96], 'trc': [415,415]})
imview(contour={'file': 'CO21PV.image', 'levels': [-25, -5, -1, 1,2,3,4,5], 'unit': 0.004})
exportfits('CO21PV.image', fitsimage='CO21PV.fits', velocity=True, optical=False, minpix=0, maxpix=0.01, dropstokes=True)
