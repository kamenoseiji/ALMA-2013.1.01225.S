prefixList =  ['uid___A002_Xa7a216_X1c4c','uid___A002_Xa7c533_X1dc7']
#uid___A002_Xa7a216_X1c4c
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0238+1636          02:38:38.930110 +16.36.59.27460 J2000   0        3116369
#  1    none J0238+166           02:38:38.930107 +16.36.59.27462 J2000   1         945952
#  2    none J0241-0815          02:41:04.798500 -08.15.20.75180 J2000   2        1092363
#  3    none J0301+0118          03:01:23.606990 +01.18.35.99630 J2000   3         965263
#  4    none NGC_1052            02:41:04.798510 -08.15.20.75170 J2000   4       10952863
#
#uid___A002_Xa7c533_X1dc7
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0339-0146          03:39:30.937790 -01.46.35.80410 J2000   0        3268104
#  1    none J0334-401           03:34:13.654489 -40.08.25.39788 J2000   1        1537494
#  2    none J0241-0815          02:41:04.798500 -08.15.20.75180 J2000   2        1690962
#  3    none J0301+0118          03:01:23.606990 +01.18.35.99630 J2000   3        1012242
#  4    none NGC_1052            02:41:04.798510 -08.15.20.75170 J2000   4       11485740
BPCals = ['J0238+1636', 'J0339-0146']
PHCals1= ['J0238+1636', 'J0334-401']
PHCals2= ['J0301+0118', 'J0301+0118']
TARGET = 'NGC_1052'
REFANT = 'DA55'
for prefix in prefixList:
    os.system('rm -rf ' + prefix + '_flagonline.txt')
    os.system('rm -rf ' + prefix + '.listobs')
    importasdm(prefix)
    listobs(vis=prefix+'.ms', scan='', spw='', verbose=True, listfile=prefix+'.listobs')
    plotants( vis=prefix + '.ms', figfile=prefix + '_plotants.png')
    browsetable(tablename = prefix +'.ms')
#
#-------- Tsys
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    flagdata(vis=prefix + '.ms', mode='manual', spw='1~24', autocorr=True, flagbackup = False)
    flagdata(vis=prefix + '.ms', mode='manual', intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*', flagbackup = False)
    flagmanager(vis=prefix+'.ms', mode='save', versionname='Apriori')
    os.system('rm -rf cal.' + prefix + '.tsys')
    gencal(vis=prefix+'.ms', caltype='tsys', caltable='cal.'+prefix+'.tsys')
    plotcal(caltable='cal.'+prefix+'.tsys', xaxis='freq', yaxis='tsys', spw='9:4~123,11:4~123,13:4~123,15:4~123', iteration='antenna', subplot=221)
#
#-------- WVR
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    os.system('rm -rf cal.' + prefix + '.WVR')
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    wvrgcal(segsource=True, caltable='cal.'+prefix+'.WVR', vis=prefix+'.ms', wvrspw=[0], tie=[BPCAL + ',' + PHCAL1+ ',' + PHCAL2+ ',' + TARGET], toffset=0, wvrflag=[], statsource=TARGET)
    plotcal('cal.'+prefix+'.WVR', xaxis='time',yaxis='phase',iteration='antenna', subplot=221)
#
#-------- Apply Tsys and WVR correction
from recipes.almahelpers import tsysspwmap
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable='cal.'+prefix+'.tsys', tsysChanTol=1)
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    for fields in list(set([BPCAL, PHCAL1, PHCAL2, TARGET])):
        applycal(vis=prefix+'.ms', field=fields, flagbackup=False, spw='17,19,21,23', interp=['linear', 'nearest'], gaintable=['cal.'+prefix+'.tsys', 'cal.'+prefix+'.WVR'], gainfield=[fields, fields], spwmap=[tsysmap, []], calwt=True)
#
#-------- Phase Cal
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    #-------- Split
    os.system('rm -rf '+prefix+'_line.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_line.ms', datacolumn='corrected', spw='17,19,21,23')
    #-------- Phase Cal for bandpass
    os.system('rm -rf P0.'+prefix)
    gaincal(vis=prefix+'_line.ms', caltable='P0.'+prefix, spw='*', field= BPCAL, scan='', selectdata=True, solint='int', refant=REFANT, calmode='p')
    os.system('rm -rf B0.'+prefix)
    bandpass(vis = prefix + '_line.ms', caltable = 'B0.'+prefix, gaintable='P0.'+prefix, spw='*:3~124', field=BPCAL, scan='', minblperant=5, minsnr=5, solint='inf', combine='scan,field', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True, spwmap=[0,1,2,3])
    #-------- Phase Cal for all
    os.system('rm -rf P1.'+prefix)
    gaincal(vis=prefix+'_line.ms', caltable='P1.'+prefix, spw='*:3~124', gaintable = ['B0.'+prefix, 'P0.'+prefix], field='', selectdata=True, solint='int', refant=REFANT, gaintype='G', combine='spw', calmode='p', minsnr=3)
    plotcal(caltable = 'P1.'+prefix, xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-180,180], iteration = 'antenna', figfile='cal_phase.png', subplot = 221)
#
#-------- Flux cal
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    os.system('rm -rf G0.'+prefix)
    setjy( vis = prefix + '_line.ms', field=TARGET, spw='0,1,2,3', standard='manual', fluxdensity=[0.4420919, 0, 0, 0], spix = [-0.6,0], reffreq = '349.798GHz', usescratch=False)
    gaincal(vis = prefix + '_line.ms', caltable = 'G0.'+prefix, spw ='*', field = '', minsnr=5.0, solint='inf', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['B0.'+prefix,'P0.'+prefix,'P1.'+prefix], spwmap=[[0,1,2,3],[0,1,2,3],[0,0,0,0]], calmode = 'a')
    plotcal(caltable = 'G0.'+prefix, xaxis = 'time', yaxis = 'amp', plotsymbol='o', plotrange = [], iteration = 'antenna', figfile='G0.'+prefix+'.png', subplot = 221)
    fluxscale(vis= prefix + '_line.ms', caltable='G0.'+prefix, fluxtable='G0.'+prefix+'.flux', reference=TARGET, transfer=BPCAL + ', ' +  PHCAL1 + ', ' + PHCAL2, refspwmap=[0,1,2,3])
    applycal(vis= prefix + '_line.ms', flagbackup=False, field='', interp=['nearest','nearest','nearest','nearest'], gainfield = '', gaintable=['B0.'+prefix, 'P0.'+prefix, 'P1.'+prefix, 'G0.'+prefix+'.flux'], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]])
#
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=0 (freq=3.44131e+11 Hz) is: 3.14133 +/- 0.00611614 (SNR = 513.613, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=1 (freq=3.43255e+11 Hz) is: 3.14797 +/- 0.00680372 (SNR = 462.684, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=2 (freq=3.55016e+11 Hz) is: 3.08471 +/- 0.0243695 (SNR = 126.581, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=3 (freq=3.57012e+11 Hz) is: 3.07054 +/- 0.00821443 (SNR = 373.798, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=0 (freq=3.44131e+11 Hz) is: 3.24987 +/- 0.00613551 (SNR = 529.682, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=1 (freq=3.43255e+11 Hz) is: 3.2571 +/- 0.0064195 (SNR = 507.376, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=2 (freq=3.55016e+11 Hz) is: 3.19992 +/- 0.0255469 (SNR = 125.257, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0238+1636 in SpW=3 (freq=3.57012e+11 Hz) is: 3.19433 +/- 0.00763009 (SNR = 418.649, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0301+0118 in SpW=0 (freq=3.44131e+11 Hz) is: 0.351933 +/- 0.00812889 (SNR = 43.2941, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0301+0118 in SpW=1 (freq=3.43255e+11 Hz) is: 0.348748 +/- 0.00486895 (SNR = 71.627, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0301+0118 in SpW=2 (freq=3.55016e+11 Hz) is: 0.399119 +/- 0.056529 (SNR = 7.06044, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Flux density for J0301+0118 in SpW=3 (freq=3.57012e+11 Hz) is: 0.336929 +/- 0.00574407 (SNR = 58.6569, N = 82)
#2019-09-20 20:51:54 INFO fluxscale	 Fitted spectrum for J0238+1636 with fitorder=1: Flux density = 3.11019 +/- 0.000691745 (freq=349.799 GHz) spidx: a_1 (spectral index) =-0.625027 +/- 0.0122622 covariance matrix for the fit:  covar(0,0)=1.9969e-06 covar(0,1)=0.000116075 covar(1,0)=0.000116075 covar(1,1)=0.0321813
#2019-09-20 20:51:54 INFO fluxscale	 Fitted spectrum for J0238+1636 with fitorder=1: Flux density = 3.22576 +/- 0.0011203 (freq=349.799 GHz) spidx: a_1 (spectral index) =-0.485416 +/- 0.0190168 covariance matrix for the fit:  covar(0,0)=1.67119e-06 covar(0,1)=8.91783e-05 covar(1,0)=8.91783e-05 covar(1,1)=0.0265662
#2019-09-20 20:51:54 INFO fluxscale	 Fitted spectrum for J0301+0118 with fitorder=1: Flux density = 0.343945 +/- 0.00303755 (freq=349.799 GHz) spidx: a_1 (spectral index) =-0.896674 +/- 0.465728 covariance matrix for the fit:  covar(0,0)=0.000103642 covar(0,1)=0.00357455 covar(1,0)=0.00357455 covar(1,1)=1.52814
#2019-09-20 20:59:02 INFO fluxscale	 Found transfer field(s):  J0339-0146 J0334-401 J0301+0118
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0339-0146 in SpW=0 (freq=3.44131e+11 Hz) is: 2.07942 +/- 0.00611091 (SNR = 340.281, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0339-0146 in SpW=1 (freq=3.43255e+11 Hz) is: 2.07018 +/- 0.00633563 (SNR = 326.753, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0339-0146 in SpW=2 (freq=3.55016e+11 Hz) is: 2.03662 +/- 0.0282273 (SNR = 72.1508, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0339-0146 in SpW=3 (freq=3.57012e+11 Hz) is: 2.02254 +/- 0.00531125 (SNR = 380.804, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0334-401 in SpW=0 (freq=3.44131e+11 Hz) is: 0.466915 +/- 0.00326824 (SNR = 142.865, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0334-401 in SpW=1 (freq=3.43255e+11 Hz) is: 0.465564 +/- 0.0032313 (SNR = 144.08, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0334-401 in SpW=2 (freq=3.55016e+11 Hz) is: 0.456899 +/- 0.0130744 (SNR = 34.9461, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0334-401 in SpW=3 (freq=3.57012e+11 Hz) is: 0.460049 +/- 0.00263239 (SNR = 174.765, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0301+0118 in SpW=0 (freq=3.44131e+11 Hz) is: 0.199016 +/- 0.0026355 (SNR = 75.5137, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0301+0118 in SpW=1 (freq=3.43255e+11 Hz) is: 0.199307 +/- 0.00349773 (SNR = 56.9819, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0301+0118 in SpW=2 (freq=3.55016e+11 Hz) is: 0.20213 +/- 0.0215317 (SNR = 9.38756, N = 84)
##2019-09-20 20:59:03 INFO fluxscale	 Flux density for J0301+0118 in SpW=3 (freq=3.57012e+11 Hz) is: 0.195983 +/- 0.00400388 (SNR = 48.9484, N = 84)
#2019-09-20 20:59:03 INFO fluxscale	 Fitted spectrum for J0339-0146 with fitorder=1: Flux density = 2.0507 +/- 0.00351484 (freq=349.798 GHz) spidx: a_1 (spectral index) =-0.665226 +/- 0.0915955 covariance matrix for the fit:  covar(0,0)=2.72502e-06 covar(0,1)=4.18475e-05 covar(1,0)=4.18475e-05 covar(1,1)=0.0412613
#2019-09-20 20:59:03 INFO fluxscale	 Fitted spectrum for J0334-401 with fitorder=1: Flux density = 0.463288 +/- 0.000597383 (freq=349.798 GHz) spidx: a_1 (spectral index) =-0.353925 +/- 0.068502 covariance matrix for the fit:  covar(0,0)=1.37626e-05 covar(0,1)=0.000102867 covar(1,0)=0.000102867 covar(1,1)=0.205937
#2019-09-20 20:59:03 INFO fluxscale	 Fitted spectrum for J0301+0118 with fitorder=1: Flux density = 0.197759 +/- 0.000402002 (freq=349.798 GHz) spidx: a_1 (spectral index) =-0.400585 +/- 0.113045 covariance matrix for the fit:  covar(0,0)=0.000117737 covar(0,1)=0.0076447 covar(1,0)=0.0076447 covar(1,1)=1.93048
#-------- Split into target source
for fileIndex in range(len(prefixList)):
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    prefix = prefixList[fileIndex]
    os.system('rm -rf ' + TARGET + prefix + '.ms'); split(vis=prefix + '_line.ms', outputvis=TARGET+prefix+'.ms', field=TARGET, datacolumn='corrected')
    os.system('rm -rf ' + BPCAL  + prefix + '.ms'); split(vis=prefix + '_line.ms', outputvis=BPCAL+prefix+'.ms', field=BPCAL, datacolumn='corrected')
    os.system('rm -rf ' + PHCAL1 + prefix + '.ms'); split(vis=prefix + '_line.ms', outputvis=PHCAL1+prefix+'.ms', field=PHCAL1, datacolumn='corrected')
    os.system('rm -rf ' + PHCAL2 + prefix + '.ms'); split(vis=prefix + '_line.ms', outputvis=PHCAL2+prefix+'.ms', field=PHCAL2, datacolumn='corrected')
#
#-------- Concat
os.system('rm -rf ' + TARGET + '_comb.ms*')
os.system('rm -rf ' + PHCals2[0] + '_comb.ms*')
comvis = []
for prefix in prefixList: comvis.append(TARGET + prefix + '.ms')
concat(vis=comvis, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis=TARGET + '_comb.ms')
comvis = []
for prefix in prefixList: comvis.append(PHCals2[0] + prefix + '.ms')
concat(vis=comvis, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis=PHCals2[0] + '_comb.ms')
plotms(vis=TARGET + '_comb.ms',spw='0:4~124', antenna='*&', xaxis='time',yaxis='amp', avgchannel='121',coloraxis='corr',iteraxis='baseline')
#-------- Imaging calibrators
#-------- Imaging J0238+1636
os.system('rm -rf SC.Gap0')
prefix = BPCals[0] + prefixList[0]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[3.11019, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
#-------- Imaging J0339-0146
os.system('rm -rf SC.Gap0')
prefix = BPCals[1] + prefixList[1]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[2.0507, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
#-------- Imaging J0334-401
os.system('rm -rf SC.Gap0')
prefix = PHCals1[1] + prefixList[1]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='8.08s,128ch', refant=REFANT, gaintype='G', smodel=[0.463288, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='linear', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
#-------- Imaging J0301+0118
os.system('rm -rf SC.Gap0')
prefix = PHCals2[0] + '_comb'
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='8.08s,128ch', refant=REFANT, gaintype='G', smodel=[0.343945, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='linear', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=False)
#-------- Imaging NGC_1052
os.system('rm -rf SC.Gap0')
"""
prefix = TARGET + '_comb'
"""
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='0:2~6;32~50;110~124, 1:50~92;118~124, 2:4~10;90~124, 3:3~124', solint='8.08s,128ch', refant=REFANT, gaintype='G', smodel=[0.4420919, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='linear', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean(vis=prefix+'.ms', datacolumn='corrected', imagename=prefix+'_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # CO J=3-2
tclean(vis=prefix+'.ms', datacolumn='corrected', imagename=prefix+'_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='344.916247GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # SiO J=8-7 v=1
tclean(vis=prefix+'.ms', datacolumn='corrected', imagename=prefix+'_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='356.734223GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # HCO+ J=4-3 v=0
tclean(vis=prefix+'.ms', datacolumn='corrected', imagename=prefix+'_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='357.921987GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # HOC+ J=4-3 
#-------- Imaging NGC_1052 Cont
prefix = TARGET + '_comb'
os.system('rm -rf ' + prefix + '_Cont*')
tclean( vis=prefix+'.ms', datacolumn='corrected', imagename=prefix + '_ContNat', spw='0:4~6;32~50;110~123, 1:50~92;118~123, 2:4~10;90~123, 3:4~123', specmode='mfs', nterms=2, niter=1000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
imview(prefix + '_ContNat.image.pbcor')
#-------- Continuum subtraction
prefix = TARGET + '_comb'
os.system('rm -rf ' + prefix + '.ms.contsub*')
uvcontsub(vis = prefix + '.ms', spw='0', fitspw='0:4~6;32~50;110~123', solint ='int', fitorder = 0)
os.system('rm -rf ' + prefix + '_CO_Nat*')
tclean(vis=prefix+'.ms.contsub', datacolumn='corrected', imagename=prefix+'_CO_Nat', spw='0', specmode='cube', start=3, nchan=122, width=1, outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # CO J=3-2
os.system('rm -rf ' + prefix + '_CO_Bin*')
tclean(vis=prefix+'.ms.contsub', datacolumn='corrected', imagename=prefix+'_CO_Bin', spw='0', specmode='cube', start=4, nchan=30, width=4, outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # CO J=3-2
os.system('rm -rf ' + prefix + '_CO_Tap*')
tclean(vis=prefix+'.ms.contsub', datacolumn='corrected', imagename=prefix+'_CO_Tap', spw='0', specmode='cube', start=40, nchan=16, width=4, outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', uvtaper='1000klambda', threshold='0.000mJy', pbcor=True) # CO J=3-2
"""
prefix = TARGET + '_comb'
"""
os.system('rm -rf ' + prefix + '_CO_forChanMap.image*')
tclean(vis=prefix+'.ms.contsub', datacolumn='corrected', imagename=prefix+'_CO_forChanMap', spw='0', specmode='cube', start='1160km/s', nchan=16, width='45km/s', outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=512, cell='0.04arcsec', weighting='natural', uvtaper='1000klambda', threshold='0.000mJy', pbcor=True) # CO J=3-2
imview(raster={'file': prefix + '_CO_forChanMap.image.pbcor', 'range': [0.000, 0.004], 'colormap': 'Rainbow 2'}, contour={'file': prefix + '_CO_forChanMap.image', 'levels': [-100, -10, 1,2,4,8,16,32,64,128,256,512,1024], 'unit': 0.001}, zoom={'blc': [96,96], 'trc': [415,415]})
#-------- Moment-0 map
# image rms = 0.36 mJy/beam
os.system('rm -rf ' + prefix + '_CO_forChanMap.image.mom*')
immoments(prefix + '_CO_forChanMap.image', moments=[0], outfile=prefix + '_CO_forChanMap.image.mom0', chans='1~12', includepix=[7.2e-4, 100])
imview(raster={'file': prefix + '_CO_forChanMap.image.mom0', 'range': [0.0, 0.4], 'colormap': 'Rainbow 2'}, contour={'file': 'NGC_1052_comb_ContNat.image.pbcor', 'levels': [1,2,4,8,16,32,64,128,256,512,1024], 'unit': 0.0005}, zoom={'blc': [96,96], 'trc': [415,415]})
imview(raster={'file': prefix + '_CO_forChanMap.image.mom0', 'range': [0.0, 0.4], 'colormap': 'Rainbow 2'}, zoom={'blc': [96,96], 'trc': [415,415]})
#-------- Moment-1 map
os.system('rm -rf NGC_1052_comb_CO_Tap.image.mom1')
immoments('NGC_1052_comb_CO_Tap.image', moments=[1], outfile='NGC_1052_comb_CO_Tap.image.mom1', chans='1~12', includepix=[1.1e-3, 100])
imview(raster={'file': 'NGC_1052_comb_CO_Tap.image.mom1', 'range': [1200, 1800], 'colormap': 'Rainbow 3'}, contour={'file': 'NGC_1052_comb_ContNat.image.pbcor', 'levels': [1,2,4,8,16,32,64,128,256,512,1024], 'unit': 0.0005}, zoom={'blc': [64,64], 'trc': [447,447]})
#-------- PV map
os.system('rm -rf ' + prefix + '_CO_forPV*')
tclean(vis=prefix+'.ms.contsub', datacolumn='corrected', imagename=prefix+'_CO_forPV', spw='0', specmode='cube', start=41, nchan=64, width=1, outframe='LSRK', veltype='radio', restfreq='345.79599GHz', niter=10000, gain=0.1, interactive=True, imsize=384, cell='0.04arcsec', weighting='natural', threshold='0.000mJy', pbcor=True) # CO J=3-2
# PV slice (192,240) - (192, 144) with 25-pix width
imview(contour={'file': 'CO32PV.image', 'levels': [-25, -5, -1, 1,2,3,4,5], 'unit': 0.004})
exportfits('CO32PV.image', fitsimage='CO32PV.fits', velocity=True, optical=False, minpix=0, maxpix=0.01, dropstokes=True)
