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
    

    def Disaggregation(IML, Sites, ProbThresh=1.e-4):
	""" 
	Disaggregation of various parameters 
	refer to OpenSHA and Ting Lin and Jack Baker
	"""
	if self.IMR == 'CyberShake': 
	    if SiteName not in Sites.keys(): 
		print 'SiteName is not in the CyberShake Sites List'
		raise ValueError 
	    IMs = Database.ExtractCyberShakeIMs(SiteName) 
	    rate = 0.0 
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
		    condProb[sr_key] = ProbIMTlgIML
		    RupProb = RupProbs[srkey] 
		    
		    if self.SourceType == 'NonPoisson': 
			Prob1 += ProbIMTlgIML * RupProb 
		    else: 
			PoE0 *= (1-RupProb)**(ProbIMTlgIML)
		if self.SourceType == 'NonPoisson':
		    PoE0 *= (1.0 - Prob1)
	    PoE.append( 1.0 - PoE0 )

