#!/usr/bin/env python
"""
Seismic Hazard Analysis 
""" 
# local
import Forecasts 
import IMRs 
import Sites 
from utils import * 


class PSHA: 
    """
    PSHA (hazard curve and map and disaggregation)
    """ 

    def __init__(self, ERF='UCERF2', SourceType='NonPoisson', IMR='CyberShake', ERFinfo=(35,5,3,1), Ti=3.0 ):
	# Initialization 
	# (especially the earthquake rutpure forecast model, and intensity measure relations with IMT)
	
	self.ERF = ERF
        self.SourceType = SourceType    # Poisson or NonPoisson 
	self.IMR = IMR 
	self.Ti = Ti 
	if self.ERF == 'UCERF2':
	    self.rup_model_ids = ERFinfo 
	    self.ERFdatabase = Forecasts.UCERF2Database(erf_id=ERFinfo[0])
	else: 
	    # other forecast models
	    pass 

	self.RupProbs = self.ERFdatabase.ExtractRupProb() 


    def HazardCurveCalc(self, IMLs, SiteInfo, SiteData, MetaPth = './metadata' ):
	"""
	Compute Hazard Curve and Plot 
	"""
        
	IMLs = IMLs.tolist() 
	SiteName = SiteInfo['SiteName'] 
	MetaPth1 = MetaPth + '/%s'%SiteName 
	for f in [MetaPth, MetaPth1,]: 
	    if not os.path.exists( f): 
		os.mkdir(f)
	
	if self.Ti != -2.0: 
	    if self.Ti == -1.0:
		filen = 'PGA'
	    else: 
		filen = 'PSA%s'%('%.3f'%self.Ti)
	else: 
	    self.filen = 'PGV'
	MetaFile = MetaPth1 + '/Meta_HazardCurve_%s_%s_%s.py'%(self.IMR, SiteName, filen)

	if not os.path.exists( MetaFile ): 
	    PoE = []
	    CS_Sites = Sites.CyberShakeSites( SiteData )
	    if self.IMR == 'CyberShake': 
		if SiteName not in CS_Sites.keys(): 
		    print 'SiteName is not in the CyberShake Sites List'
		    raise ValueError 
		CS = IMRs.CyberShakeIMs(self.rup_model_ids, Ti=self.Ti)
		IMs = CS.ExtractIMs(SiteName) 
		for iml in xrange( len(IMLs) ):
		    IML = IMLs[iml]
		    PoE0 = 1.0 
		    for sid in IMs.keys():
			skey = '%s'%sid
			if self.SourceType == 'NonPoisson':
			    Prob1 = 0.0
			for rid in IMs[skey].keys(): 
			    srkey = '%s,%s'%(sid,rid) 
			    rkey = '%s'%rid 
			    NL = len(IMs[skey][rkey])
			    NIML = len((np.array(IMs[skey][rkey])>=IML).nonzero()[0])
			    ProbIMTlgIML = NIML*1.0/NL 
			    RupProb = self.RupProbs[srkey] 
			    if self.SourceType == 'NonPoisson': 
				Prob1 += ProbIMTlgIML * RupProb 
			    else: 
				PoE0 *= (1-RupProb)**(ProbIMTlgIML)
			if self.SourceType == 'NonPoisson':
			    PoE0 *= (1.0 - Prob1)
		    PoE.append( 1.0 - PoE0 )

	    if self.IMR in ['CB08','BA08','CY08','AS08']:
		
		# Site info 
		SiteFlag = 0
		if SiteInfo['SiteName'] in CS_Sites.keys(): 
		    print 'Given site name is in CyberShake site list'
		    SiteName = SiteInfo['SiteName']
		    SiteGeom = float(CS_Sites[SiteName]['lon']), float(CS_Sites[SiteName]['lat']),0
		    Vs30 = float(CS_Sites[SiteName]['Vs30_Wills_2006'])
		    Z10 = float(CS_Sites[SiteName]['Z1.0']) * 1000   # km to m
		    Z25 = float(CS_Sites[SiteName]['Z2.5'])
		    SiteFlag = 1
		else: 
		    try: 
			SiteGeom = SiteInfo['lon'], SiteInfo['lat'],0
		    except: 
			print 'You have to specify the loation of your site'
			raise ValueError 

		    for SiteName in CS_Sites.keys(): 
			if SiteGeom[0] == float(Sites[SiteName]['lon']) and SiteGeom[1] == float(Sites[SiteName]['lat']): 
			    print 'The site is a CyberShake site'
			    Vs30 = float(CS_Sites[SiteName]['Vs30_Wills_2006'])
			    Z10 = float(CS_Sites[SiteName]['Z1.0']) * 1000   # km to m
			    Z25 = float(CS_Sites[SiteName]['Z2.5'])
			    SiteFlag = 1
			    break
		if SiteFlag == 0:
		    try: 
			Vs30 = SiteInfo['Vs30'] 
			Z10 = None
			Z25 = None 
		    except: 
			print 'When using other IMR but CyberShake, you should give at least the Vs30 for the site!'
			raise ValueError 
		
		# Source Info (fault geometry for distance calculation) (need for NGA model)
		NGA08IM = IMRs.NGA08IMs(self.IMR, self.Ti)
		SourceInfo0 = self.ERFdatabase.ExtractSource0() 

		condProb = {}
		Nsources = len(SourceInfo0.keys())
		itime = 0
		for skey in SourceInfo0.keys():
		    Nrups = SourceInfo0[skey]
		    for rid in xrange( Nrups ):
			rkey = '%s'%rid
			SourceInfo = self.ERFdatabase.ExtractRuptures(skey,rkey)    # takes time!
			mu, std = NGA08IM.ExtractIMs(SourceInfo, SiteGeom, Vs30, Z25, Z10)
			for iml in xrange( len(IMLs) ):
			    srp_key = '%s,%s,%s'%(skey,rkey,iml)
			    IML = IMLs[iml]
			    UpperLimit = mu+4*std
			    if np.log(IML) > UpperLimit: 
				ProbIMTlgIML = 0.0
			    else: 
				ProbIMTlgIML = GaussianCumulative(np.log(IML),UpperLimit, mu, std )
			    condProb[srp_key] = ProbIMTlgIML 
		    sys.stdout.write('\rIM extraction finished %.2f%%'%(itime*100.0/Nsources))
		    sys.stdout.flush()
		    sleep(0.1)
		    itime += 1
	        sys.stdout.write('\n')

		PoE = []
		for iml in xrange( len(IMLs) ):
		    IML = IMLs[iml]
		    PoE0 = 1.0 
		    for skey in SourceInfo0.keys():
			if self.SourceType == 'NonPoisson':
			    Prob1 = 0.0
			Nrups = SourceInfo0[skey]
			for rid in xrange( Nrups ):
			    rkey = '%s'%rid
			    srkey = '%s,%s'%(skey,rkey) 
			    srp_key = '%s,%s'%(skey,rkey,iml) 
			    RupProb = self.RupProbs[srkey] 
			    ProbIMTlgIML = condProb[srp_key]   # for a given source and rupture 
			    if self.SourceType == 'NonPoisson':
				Prob1 += ProbIMTlgIML * RupProb 
			    else: 
				PoE0 *= (1-RupProb) ** ProbIMTlgIML 
			if self.SourceType == 'NonPoisson':	
			    PoE0 *= (1.0 - Prob1)
		    PoE.append( 1.0 - PoE0 )
                
	    header = '#metafile for site %s at %s sec, IMR=%s\n'%(SiteName, '%.3f'%self.Ti, self.IMR)

	    # save to MetaFile 
	    meta = dict( IMLs = IMLs, 
		         PoE = PoE,)
	    save(MetaFile, meta, header=header)

        else: 
	    meta = load( MetaFile )
	    PoE = meta.PoE 
	    IMLs = meta.IMLs 

	return PoE, IMLs 
    

    def Disaggregation(self, SourceDisagg, RuptureDisagg, IML, CyberShakeDatabase, SiteName, SiteData=None,PoE=None, IMLs=None, DisaggMeta = './DisaggMeta_CurrentStation.txt'):
	""" 
	Disaggregation of CyberShake for one site
	refer to OpenSHA and Ting Lin and Jack Baker
	"""
	stanam = SiteName

	SourceDisagg[stanam] = {}
	RuptureDisagg[stanam] = {}
	
	SourceDisagg0 = {} 
	RuptureDisagg0 = {}
	if not os.path.exists( DisaggMeta ):

	    if self.IMR == 'CyberShake': 
		if SiteData != None: 
		    CS_Sites = Sites.CyberShakeSites( SiteData )
		    if SiteName not in CS_Sites.keys(): 
			print 'SiteName is not in the CyberShake Sites List'
			raise ValueError 

		# read all IMs (both for Disagg and Hazard Curves)
		print 'IM extraction begins'
		IMs = CyberShakeDatabase.ExtractIMs(SiteName)
		print 'IM extraction ends'

		if IML != 0: 
		    
		    fid = open(DisaggMeta,'w')
		    fid.write('#SourceID RuptureID SourceDisagg RuptureDisagg\n')
		    
		    # disaggregate directly 
		    TotalRate = 0.0
		    for sid in IMs.keys():
			skey = '%s'%sid
			SourceRate = 0.0
			for rid in IMs[skey].keys(): 
			    srkey = '%s,%s'%(sid,rid) 
			    rkey = '%s'%rid 
			    
			    NL = len(IMs[skey][rkey])
			    NIML = len((np.array(IMs[skey][rkey])>=IML).nonzero()[0])
			    ProbIMTlgIML = NIML*1.0/NL 
			    RupProb = self.RupProbs[srkey] 
			    rate = -np.log( 1- RupProb ) * ProbIMTlgIML   # rate for each rupture 
			    #print 'Rate: ', rate
			    SourceRate += rate 
			    TotalRate += rate 
			    RuptureDisagg0[srkey] = rate
			SourceDisagg0[skey] = SourceRate
		    #print 'Total Rate: ', TotalRate
		    #raw_input() 

		    for skey in IMs.keys(): 
			if TotalRate == 0.0: 
			    SourceDisagg0[skey] = 0.0 
			else: 
			    SourceDisagg0[skey] = SourceDisagg0[skey] / TotalRate
			for rkey in IMs[skey].keys(): 
			    srkey = '%s,%s'%(skey,rkey)
			    if TotalRate == 0.0:
				RuptureDisagg0[srkey] = 0.0 
			    else: 
				RuptureDisagg0[srkey] /= TotalRate 
			    fid.write('%s %s %s %s\n'%(skey,rkey,SourceDisagg0[skey],RuptureDisagg0[srkey])) 
		    fid.close()

		else: 
		    if PoE != None: 
			if IMLs != None: 
			    pass
			else: 
			    print 'You should specify IMLs list to compute Hazard Curve' 
			    raise ValueError
		    else: 
			print 'You should specify PoE when IML = 0' 
			raise ValueError 
	else: 
	    # read from file :
	    # attention to those nan (total rate == 0)
	    data = np.loadtxt( DisaggMeta, skiprows=1 ) 
	    for irow in xrange( len(data) ): 
		sid,rid,tmp1,tmp2 = data[irow]
		if np.isnan(tmp1):
		    # when total rate is 0.0, this will give you nan /TotalRate
		    # the reason is that your IML is too high for specified site (all source will not make to IML)  
		    # In this situation, the disaggregation will be 0.0
		    tmp1 = 0.0 
		    tmp2 = 0.0
		skey = '%s'%(int(sid)) 
		rkey = '%s'%(int(rid)) 
		srkey = '%s,%s'%(skey,rkey)
		SourceDisagg0[skey] = tmp1 
		RuptureDisagg0[srkey] = tmp2
	
	SourceDisagg[stanam] = SourceDisagg0 
	RuptureDisagg[stanam] = RuptureDisagg0 
       
        return SourceDisagg, RuptureDisagg 

