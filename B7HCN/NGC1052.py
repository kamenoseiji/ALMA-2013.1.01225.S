prefixList =  ['uid___A002_Xa830fc_X16e2','uid___A002_Xa830fc_X1a2f']
#uid___A002_Xa830fc_X16e2
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0224+0659          02:24:28.428200 +06.59.23.34160 J2000   0        2412324
#  1    none J2357-5311          23:57:53.266077 -53.11.13.68929 J2000   1        1135044
#  2    none J0241-0815          02:41:04.798500 -08.15.20.75180 J2000   2        1006848
#  3    none J0219+0120          02:19:07.024510 +01.20.59.86610 J2000   3         498024
#  4    none NGC_1052            02:41:04.798510 -08.15.20.75170 J2000   4        5935140
#
#uid___A002_Xa830fc_X1a2f
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0423-0120          04:23:15.800730 -01.20.33.06550 J2000   0        2412396
#  1    none J0238+166           02:38:38.930107 +16.36.59.27462 J2000   1        1135044
#  2    none J0241-0815          02:41:04.798500 -08.15.20.75180 J2000   2        1087416
#  3    none J0219+0120          02:19:07.024510 +01.20.59.86610 J2000   3         996300
#  4    none NGC_1052            02:41:04.798510 -08.15.20.75170 J2000   4       11870568
BPCals = ['J0224+0659', 'J0423-0120']
PHCals1= ['J2357-5311', 'J0238+166']
PHCals2= ['J0219+0120', 'J0219+0120']
TARGET = 'NGC_1052'
REFANT = 'DA55'
"""
for prefix in prefixList:
    os.system('rm -rf ' + prefix + '_flagonline.txt')
    importasdm(prefix)
    listobs(vis=prefix+'.ms', scan='', spw='', verbose=True, listfile=prefix+'.listobs')
    plotants( vis=prefix + '.ms', figfile=prefix + '_plotants.png')
    browsetable(tablename = prefix +'.ms')
#
from recipes.almahelpers import tsysspwmap
for fileIndex in range(len(prefixList)):
    prefix = prefixList[fileIndex]
    BPCAL = BPCals[fileIndex]
    PHCAL1 = PHCals1[fileIndex]
    PHCAL2 = PHCals2[fileIndex]
    #-------- Flagging a priori
    #flagdata(vis=prefix + '.ms', mode='summary', name='after')
    flagdata(vis=prefix + '.ms', mode='manual', spw='1~24', autocorr=True, flagbackup = False)
    flagdata(vis=prefix + '.ms', mode='manual', intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*', flagbackup = False)
    flagmanager(vis=prefix+'.ms', mode='save', versionname='Apriori')
    #-------- Tsys and WVR
    os.system('rm -rf cal.' + prefix + '.WVR')
    os.system('rm -rf cal.' + prefix + '.tsys')
    gencal(vis=prefix+'.ms', caltype='tsys', caltable='cal.'+prefix+'.tsys')
    if fileIndex == 0: au.copyTsys('cal.'+prefix+'.tsys', fromAntenna='DA60', toAntenna='DA62', spw=[9, 11, 13, 15], scan=[3, 6, 10, 12, 17, 19], scaleFactor=1.0, unflag=True)
    if fileIndex == 1: au.copyTsys('cal.'+prefix+'.tsys', fromAntenna='DA60', toAntenna='DA62', spw=[9, 11, 13, 15], scan=[3, 6, 9, 11, 16, 18, 23, 25, 30, 32], scaleFactor=1.0, unflag=True)
    plotcal(caltable='cal.'+prefix+'.tsys', xaxis='freq', yaxis='tsys', spw='9:4~123, 11:4~123, 13:4~123, 15:4~123', iteration='antenna', subplot=221)
    os.system('rm -rf cal.' + prefix + '.WVR')
    wvrgcal(segsource=True, caltable='cal.'+prefix+'.WVR', vis=prefix+'.ms', wvrspw=[0], tie=[BPCAL + ',' + PHCAL1+ ',' + PHCAL2+ ',' + TARGET], toffset=0, wvrflag=[], statsource=TARGET)
    plotcal('cal.'+prefix+'.WVR', xaxis='time',yaxis='phase',iteration='antenna', subplot=221)
    applycal(vis=prefix+'.ms', flagbackup=False, spw='17,19,21,23', interp='linear', gaintable='cal.'+prefix+'.WVR')
    tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable='cal.'+prefix+'.tsys', tsysChanTol=1)
    for fields in list(set([BPCAL, PHCAL1, PHCAL2, TARGET])):
        applycal(vis=prefix+'.ms', field=fields, flagbackup=False, spw='17,19,21,23', interp='linear', gaintable='cal.'+prefix+'.tsys', gainfield=fields, spwmap=tsysmap, calwt=True)
    #-------- Split
    os.system('rm -rf '+prefix+'_line.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_line.ms', datacolumn='corrected', spw='17,19,21,23')
    os.system('rm -rf '+prefix+'_chav.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_chav.ms', datacolumn='corrected', spw='18,20,22,24')
#
#-------- Concat
os.system('rm -rf Comb_line.ms*')
comLine, comChav = [], []
for prefix in prefixList:
    comLine.append(prefix + '_line.ms')
    comChav.append(prefix + '_chav.ms')
#
concat(vis=comLine, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis= 'Comb_line.ms')
concat(vis=comChav, freqtol='0.5MHz', dirtol='0.1arcsec', concatvis= 'Comb_chav.ms')
listobs(vis='Comb_line.ms', scan='', spw='', verbose=True, listfile='Comb.listobs')
#-------- Phase Cal for bandpass
os.system('rm -rf P0')
gaincal(vis='Comb_chav.ms', caltable='P0', spw='*', field= BPCals[0] + ',' + BPCals[1] + ',' + PHCals1[1], scan='', selectdata=True, solint='int', refant=REFANT, calmode='p')
plotcal(caltable = 'P0', xaxis = 'time', yaxis = 'phase', plotsymbol='o', plotrange = [], iteration = 'antenna', figfile='P0.png', subplot = 221)
os.system('rm -rf B0')
bandpass(vis = 'Comb_line.ms', caltable = 'B0', gaintable='P0', spw='*:3~124', field=BPCals[0] + ',' + BPCals[1] + ',' + PHCals1[1], scan='4,24,27', minblperant=5, minsnr=5, solint='inf', combine='obs,scan,field', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True, spwmap=[0,1,2,3])
plotbandpass(caltable = 'B0', xaxis='freq', yaxis='amp', plotrange = [0,0,0,1.2], figfile='B0.png')
plotbandpass(caltable = 'B0', xaxis='freq', yaxis='phase',plotrange = [0,0,-180.0,180.0], figfile='B0.png')
#-------- Phase Cal to align SPW 
os.system('rm -rf P1')
gaincal(vis='Comb_line.ms', caltable='P1', spw='*:3~124', gaintable = ['B0','P0'], field=BPCAL[2], selectdata=True, solint='inf', refant=REFANT, gaintype='G', calmode='p', minsnr=7, spwmap=[[0,1,2,3],[0,0,0,0]])
plotcal(caltable = 'P1', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-180,180], iteration = 'antenna', figfile='P1.png', subplot = 221)
#-------- Phase Cal for all
os.system('rm -rf P2')
gaincal(vis='Comb_line.ms', caltable='P2', spw='*:3~124', gaintable = ['B0','P1'], field='0,1,2,3,4,5,6', selectdata=True, solint='4.04s', refant=REFANT, gaintype='G', combine='spw', calmode='p', minsnr=3)
plotcal(caltable = 'P2', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-180,180], iteration = 'antenna', figfile='P2.png', subplot = 221)
#-------- Flux cal
os.system('rm -rf G0')
setjy( vis = 'Comb_line.ms', field=TARGET, spw='0,1,2,3', standard='manual', fluxdensity=[0.4420919, 0, 0, 0], spix = [-0.6,0], reffreq = '349.798GHz', usescratch=False)
gaincal(vis = 'Comb_line.ms', caltable = 'G0', spw ='', field = '', minsnr=5.0, solint='inf', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['B0','P1','P2'], spwmap=[[0,1,2,3],[0,1,2,3],[0,0,0,0]], calmode = 'a')
plotcal(caltable = 'G0', xaxis = 'time', yaxis = 'amp', plotsymbol='o', plotrange = [], iteration = 'antenna', figfile='G0.png', subplot = 221)
fluxscale(vis= 'Comb_line.ms', caltable='G0', fluxtable='G0.flux', reference=TARGET, transfer='0,1,3,5,6', refspwmap=[0,1,2,3])
applycal(vis= 'Comb_line.ms', flagbackup=False, field='', interp=['nearest','nearest','nearest','nearest'], gainfield = '', gaintable=['B0', 'P1', 'P2', 'G0.flux'], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]])
plotms(vis='Comb_line.ms',spw='0:4~124', antenna='*&', xaxis='time',yaxis='amp', avgchannel='121',coloraxis='corr',iteraxis='baseline')
#BPCals = ['J0224+0659', 'J0423-0120']
#PHCals1= ['J2357-5311', 'J0238+166']
#PHCals2= ['J0219+0120', 'J0219+0120']
#TARGET = 'NGC_1052'
BPCals[0]  J0224+0659 0.873914 Jy
PHCals1[0] J2357-5311 0.535532 Jy
PHCals2[0] J0219+0120 0.174406 Jy
BPCals[1]  J0423-0120 1.22278  Jy
PHCals1[1] J0238+166  3.37334  Jy
#-------- Split into target source
os.system('rm -rf ' + BPCals[0] + '.ms'); split(vis='Comb_line.ms', outputvis=BPCals[0]+'.ms', field=BPCals[0], datacolumn='corrected')
os.system('rm -rf ' + BPCals[1] + '.ms'); split(vis='Comb_line.ms', outputvis=BPCals[1]+'.ms', field=BPCals[1], datacolumn='corrected')
os.system('rm -rf ' + PHCals1[0] + '.ms'); split(vis='Comb_line.ms', outputvis=PHCals1[0]+'.ms', field=PHCals1[0], datacolumn='corrected')
os.system('rm -rf ' + PHCals1[1] + '.ms'); split(vis='Comb_line.ms', outputvis=PHCals1[1]+'.ms', field=PHCals1[1], datacolumn='corrected')
os.system('rm -rf ' + PHCals2[0] + '.ms'); split(vis='Comb_line.ms', outputvis=PHCals2[0]+'.ms', field=PHCals2[0], datacolumn='corrected')
os.system('rm -rf ' + TARGET + '.ms'); split(vis='Comb_line.ms', outputvis=TARGET+'.ms', field=TARGET, datacolumn='corrected')
#-------- Imaging NGC_1052 Cont
"""
prefix = TARGET
os.system('rm -rf ' + prefix + '_Cont*')
tclean(vis=prefix+'.ms', datacolumn='corrected', imagename=prefix + '_Cont', spw='0:4~32;86~123, 1:4~89, 2:4~35;60~115, 3:30~47;100~123', specmode='mfs', nterms=2, niter=1000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
"""
#-------- Imaging NGC 1052
os.system('rm -rf SC.Gap0')
prefix = TARGET
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='0:4~32;86~123, 1:4~89, 2:4~35;60~115, 3:30~47;100~123', solint='4.04s,128ch', refant=REFANT, gaintype='G', smodel=[0.4420919, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True, restfreq='354.505473GHz') # HCN v=0 J=4-3
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True, restfreq='353.622753GHz') # H26alpha
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True, restfreq='342.88285GHz') # CS J=7-6
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=512, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True, restfreq='340.71416GHz') # SO J=7_8-6_7
#-------- Imaging J0238+1660
os.system('rm -rf SC.Gap0')
prefix = PHCals1[1]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='int,128ch', refant=REFANT, gaintype='G', smodel=[3.37334, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
#-------- Imaging J0423-0120
os.system('rm -rf SC.Gap0')
prefix = BPCals[1]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='4.04s,128ch', refant=REFANT, gaintype='G', smodel=[1.222786, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
#-------- Imaging J0219+0120
os.system('rm -rf SC.Gap0')
prefix = PHCals2[0]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='8.08s,128ch', refant=REFANT, gaintype='G', smodel=[0.174406, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
#-------- Imaging J2357-5311
os.system('rm -rf SC.Gap0')
prefix = PHCals1[0]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='4.04s,128ch', refant=REFANT, gaintype='G', smodel=[0.535532, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
#-------- Imaging J0224+0659
os.system('rm -rf SC.Gap0')
prefix = BPCals[0]
gaincal(vis=prefix+'.ms', caltable='SC.Gap0', spw='*:3~124', solint='4.04s,128ch', refant=REFANT, gaintype='G', smodel=[0.873914, 0.0,  0.0, 0], calmode='ap', minblperant=5, minsnr=3)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'phase', plotsymbol='.', plotrange = [0,0,-20,20], iteration = 'antenna', figfile='SC.Gap0.phase.png', subplot = 221)
plotcal(caltable = 'SC.Gap0', xaxis = 'time', yaxis = 'amp', plotsymbol='.', plotrange = [0,0,0.0,2.0], iteration = 'antenna', figfile='SC.Gap0.amp.png', subplot = 221)
applycal(vis=prefix+'.ms', interp='nearest', gaintable='SC.Gap0', flagbackup=False)
os.system('rm -rf ' + prefix + '_SPW*')
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW0', spw='0', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW1', spw='1', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW2', spw='2', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)
tclean( vis=prefix+'.ms', imagename=prefix + '_SPW3', spw='3', specmode='cube', start=3, nchan=122, width=1, niter=10000, threshold='0.000mJy', imsize=256, gain=0.1, cell='0.04arcsec', weighting='natural', interactive=True, pbcor=True)

