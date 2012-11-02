#!/usr/bin/env python
"""
Utilities for SHA python version 
""" 
import os, sys
import numpy as np 

import MySQLdb as mdb

from pynga.utils import * 
from pynga import * 

# intensity measure types 
IMTdict = {'-1.0':'PGA', '-2.0':'PGV', 'T':'PSA'}


def time_integral( ft, a, b, h, met = 'Trapezoid' ):
    """
    Compute numerical integral
    Input:
	ft: input time serise
	a: lower integration limit
	b: higher integration limit
	h: time interval
	met: method to do the numerical integration
    Output:
        integral
    """
    #nt1 = int( a/h + 1.0 ); nt2 = int( b/h + 1.0 )
    #f = ft[nt1:nt2]
    f = ft
    sum = 0.0 

    if met == 'Trapezoid':

	n = int( ( b - a ) / h +1)
	sum = 0.5 * ( f[0] + f[n-1] ) 
	
	for i in xrange( 1, n-1 ):
	    sum = sum + f[i]
	
	sum = sum * h

    elif met == 'Simpson':
        
	n = int( (b-a)/h+1 )

	sum_odd = f[1]
	sum_even = 0.0
	sum = f[0] + f[n-1]
	
	for i in xrange( 1, int( n/2 ) ): 
	    sum_even = sum_even + f[2*i]
	    sum_odd = sum_odd + f[2*i+1]
        
	sum_odd = 4.*sum_odd
	sum_even = 2.*sum_even

	sum = h/3 * ( sum + sum_odd + sum_even )
    
    return sum


def GaussianDistribution(x,mu=0,sigma=1):
    A = x-mu
    B = 1./(sigma*np.sqrt(2.*np.pi))
    C = np.exp( -0.5*(A/sigma)**2. )
    return B*C

def GaussianCumulative(LowerLimit, UpperLimit, mu,std, N = 100, met='Simpson'): 
    h = (UpperLimit-LowerLimit) / (N-1)
    x = np.arange( N ) * h + LowerLimit
    ft = GaussianDistribution(x, mu=mu, sigma=std) 
    integral = time_integral( ft, LowerLimit, UpperLimit, h, met=met )
    return integral 

def CyberShakeSites(sitedata, VelModel='CVMS'):
    """
    CyberShake is site-based model
    All CyberShake Sites information (Vs30, Z2.5, Z1.0 etc.)
    sitedata give the full path of the site file
    """
    sites = {}
    if VelModel == 'CVMS':
	# required header of the sitedata file
	keys = 'lat','lon', 'Vs30_Wills_2006', 'Vs30_Wald_2007', 'Z2.5', 'Z1.0'  # Z2.5,Z1.0 all in km !!
	lines = open( sitedata, 'r' ).readlines()
	for i in range( 1, len(lines) ):
	    spl = lines[i].strip().split()
	    name = spl[0]
	    sites[name] = {}
	    for ikey in xrange( len(keys) ):
		sites[name][keys[ikey]] = spl[ikey+1]
    else: 
	# For other CVM models (ask kevin for generating site condition files for CVMH model)
	#...
	print 'Other site conditions'
    return sites




# CyberShake Database (extract PSA, and probability of ruptures)
class CyberShakeDatabase:
    """
    CyberShake database 
    Site, Sources for SHA analysis
    """
    def __init__(self,rup_model_ids,Ti=3.0,Debug=False):
	"""
	Initiation of CyberShake database
	wrk: where to use CyberShake database
	rup_model_ids: CyberShake study
	"""
	# CyberShake database
	hostname = 'focal.usc.edu'  # host                            
	username = 'cybershk_ro'    # usrname                         
	password = 'CyberShake2007' # password                        
	database = 'CyberShake'     # Database name 

	# connect to the database
	try:
	    self.db = mdb.connect( host = hostname, user = username, \
		    passwd = password, db = database )
	except:
	    self.db = None
	
	self.erf_id, self.sgt_id, self.rup_scenario_id, self.vel_id = rup_model_ids
	self.Ti = Ti 
        
	#debug flat 
	self.Debug = Debug



    def ExtractRupProb(self): 
	
	cursor = self.db.cursor()
	query = 'Select Source_ID, Rupture_ID, Prob from Ruptures where ERF_ID=%s'%self.erf_id 
	cursor.execute( query ) 
	rows = cursor.fetchall()
	ERFprob = {}
	for irow in xrange( len( rows ) ):
	    sid = rows[irow][0]
	    rid = rows[irow][1]
	    srkey = '%s,%s'%(sid,rid)
	    prob = rows[irow][2]
	    ERFprob[srkey] = prob 
        cursor.close() 
	return ERFprob



    def ExtractSources(self): 
	
	SourceInfo = {}
	
	cursor = self.db.cursor()
	query = 'Select Source_ID from Ruptures where ERF_ID=%s and Rupture_ID=0'%self.erf_id 
	cursor.execute( query ) 
	rows = cursor.fetchall()
	for irow in xrange( len(rows) ): 
	    sid = rows[irow][0] 
	    skey = '%s'%sid
	    SourceInfo[skey] = {} 
	    
	    tmp = {} 
	    # Ruptures
	    query = "select * from %s where %s = %s and %s = %s"%\
		    ('Ruptures','ERF_ID',self.erf_id,'Source_ID',sid)
	    cursor.execute(query)
	    row_rup = cursor.fetchall() 

	    tmp = {}
	    SourceName = row_rup[0][3]   # section name
	    SourceInfo[skey]['Name'] = SourceName 

	    Mws = []; rids = []
	    for irup in xrange( len( row_rup ) ):
		rid = row_rup[irup][2]
		Mw = row_rup[irup][5] 
	        tmp['%s'%rid] = Mw 
	    
	    SourceInfo[skey]['Ruptures'] = tmp 
	    
	    # Fault Trace
	    query = "select * from %s where %s = %s and %s = %s and Rupture_ID=0"%('Points','ERF_ID',self.erf_id,'Source_ID',sid)
	    cursor.execute( query )       # run query
	    row_point = cursor.fetchall()
	    nvar = len(row_point)   
	    
	    Ztor = row_point[0][6]
	    Zbom = row_point[-1][6] 

	    rake = 0.0; dip = 0.0 
	    lon = []; lat = []
	    for i in xrange( len(row_point) ):
		rake = rake + row_point[i][7]
		dip = dip +  row_point[i][8]
		if row_point[i][6]==Ztor:
		    lon.append( row_point[i][5] )
		    lat.append( row_point[i][4] )

	    # compute averaged rake and dip angle
	    rake = rake / (i+1)
	    dip = dip / (i+1)
	    
	    SourceInfo[skey]['FaultInfo'] = [lon,lat,dip,Ztor,Zbom,rake] 

        return SourceInfo 


    def ExtractCyberShakeIMs(self,CS_Short_Name):
	
	cursor = self.db.cursor()
	
	# get site id
	query = "select CS_Site_ID from CyberShake_Sites where CS_Short_Name='%s'"%(CS_Short_Name)
	cursor.execute( query )
	Site_ID = cursor.fetchone()[0]

	# get run id
	query = "select Run_ID from CyberShake_Runs where Site_ID = %s and Rup_Var_Scenario_ID=%s and Velocity_Model_ID = %s and ERF_ID= %s and SGT_Variation_ID=%s and Status = 'Verified'"%(Site_ID,self.rup_scenario_id,self.vel_id,self.erf_id,self.sgt_id)
	cursor.execute( query )
	try: 
	    Run_ID = cursor.fetchone()[0]
	except: 
	    print 'There is no verified runs for your parameter list'
	    raise ValueError 
	
	# get psa at specified period
	query = "select * from %s"%('IM_Types')
	cursor.execute( query )       # run query
	row_imtype = cursor.fetchall()  
	for ir in xrange( len(row_imtype) ):
	    if abs( self.Ti - row_imtype[ir][2] ) < 1.e-4:
		IM_Type_ID = row_imtype[ir][0] 
		break

	query = "select Source_ID from CyberShake_Site_Ruptures where CS_Site_ID=%s and ERF_ID = %s and Rupture_ID = 0"%(Site_ID,self.erf_id)
	cursor.execute( query )
	Source_IDs = cursor.fetchall()
	
	if self.Debug:
	    # simulated sources and ruptures
	    query = "select Source_ID from CyberShake_Site_Ruptures where CS_Site_ID=%s and ERF_ID = %s and Rupture_ID = 0"%(Site_ID,self.erf_id)
	    cursor.execute( query )
	    Source_IDs = cursor.fetchall()
	    Nsrc = len(Source_IDs)

	    query = "select Rupture_ID from CyberShake_Site_Ruptures where CS_Site_ID=%s and ERF_ID = %s"%(Site_ID,self.erf_id)
	    cursor.execute( query )
	    Rupture_IDs = cursor.fetchall()
	    Nrup = len(Rupture_IDs)
	    print 'Source: %s, Rupture: %s'%(Nsrc,Nrup)

	    # all available PeakAmplitudes
	    query = "select * from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s "%(Run_ID,IM_Type_ID)
	    cursor.execute( query )
	    Npsa0 = len( cursor.fetchall() )
	    
	    Nvar = 0; Npsa = 0; Ndiff = 0

	if self.Ti != -2.0: 
	    factor = 1./980.  # convert PGA and SA to (g)
	else: 
	    factor = 1.0    # keep PGV at cm/s

        IMs = {} 
	for isrc in xrange( len(Source_IDs) ):
	    sid = Source_IDs[isrc][0]
	    skey = '%s'%sid 
	    IMs[skey] = {} 

	    query = "select Rupture_ID from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s and Source_ID=%s and Rup_Var_ID=0"%(Run_ID,IM_Type_ID,sid )
	    cursor.execute( query )
	    Rupture_IDs = cursor.fetchall()
	    for irup in xrange( len(Rupture_IDs) ):
		rid = Rupture_IDs[irup][0]
		rkey = '%s'%rid 
		IMs[skey][rkey] = []
		query = "select IM_Value from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s and Source_ID=%s and Rupture_ID=%s"%(Run_ID,IM_Type_ID,sid,rid)
		cursor.execute( query )
		IMvalues = cursor.fetchall() 
		for im in xrange( len(IMvalues) ): 
		    IMs[skey][rkey].append( IMvalues[im][0]*factor)     # unit conversion
	        
	    if self.Debug:
		query = "select IM_Value from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s and Source_ID=%s"%(Run_ID,IM_Type_ID,sid)
		cursor.execute( query )
		Ntmp_PA = len( cursor.fetchall() )
		Npsa = Npsa + Ntmp_PA
		query = "select * from Rupture_Variations where ERF_ID=%s and Rup_Var_Scenario_ID=%s and Source_ID = %s"%(self.erf_id, self.rup_scenario_id,Source_IDs[isrc][0])
		cursor.execute( query )
		Ntmp_RV = len( cursor.fetchall() )
		if Ntmp_PA != Ntmp_RV: 
		    # test rupture number and rupture variation number that cause the difference 
		    query = "select Rupture_ID from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s and Source_ID=%s and Rup_Var_ID=0"%(Run_ID,IM_Type_ID,Source_IDs[isrc][0])
		    cursor.execute( query )
		    Nrup_PA = len(cursor.fetchall())

		    query = "select Rupture_ID from Rupture_Variations where ERF_ID=%s and Rup_Var_Scenario_ID=%s and Source_ID = %s and Rup_Var_ID=0"%(self.erf_id, self.rup_scenario_id,Source_IDs[isrc][0])
		    cursor.execute( query )
		    Nrup_RV = len( cursor.fetchall() ) 
		    
		    query = "select Rup_Var_ID from PeakAmplitudes where Run_ID=%s and IM_Type_ID=%s and Source_ID=%s and Rupture_ID=0"%(Run_ID,IM_Type_ID,Source_IDs[isrc][0])
		    cursor.execute( query )
		    Nvar_PA = len(cursor.fetchall())

		    query = "select Rup_Var_ID from Rupture_Variations where ERF_ID=%s and Rup_Var_Scenario_ID=%s and Source_ID = %s and Rupture_ID=0"%(self.erf_id, self.rup_scenario_id,Source_IDs[isrc][0])
		    cursor.execute( query )
		    Nvar_RV = len( cursor.fetchall() ) 
		    print 'Source %s'%Source_IDs[isrc][0]
		    print 'PSA values: %s, RV values: %s'%(Ntmp_PA,Ntmp_RV)
		    print 'Rupture PA: %s, Rupture RV: %s'%(Nrup_PA,Nrup_RV)
		    print 'Rvar PA: %s, Rvar RV: %s'%(Nvar_PA, Nvar_RV)
        cursor.close() 

	return IMs


class PSHA: 
    """
    PSHA (hazard curve and map and disaggregation)
    """ 
    def __init__(self,SiteInfo,ERFinfo,Ti, IMR='CyberShake',ERF='UCERF2', MetaPth='./'):

	self.IMR = IMR 
	self.ERF = ERF
	if self.ERF == 'UCERF2':
	    self.rup_model_ids = ERFinfo 
	self.Ti = Ti 

	if self.IMR == 'CyberShake': 
	    self.SiteName = SiteInfo['SiteName'] 
	else: 
	    # you need to read rupture surface in formation to compute GMPEs values and use Gaussian (refer the following block of code)
	    self.SiteInfo = SiteInfo   # dictionary with SiteName, lon, lat, Vs30 [Z1.0, Z2.5] For GMPEs 
	    self.SiteName = SiteInfo['SiteName'] 

	if self.Ti != -2.0: 
	    if self.Ti == -1.0:
		self.xlab = 'PGA (g)'
		self.filen = 'PGA'
	    else: 
		self.xlab = '%s sec PSA (g)'%('%.3f'%self.Ti)
		self.filen = 'PSA%s'%('%.3f'%self.Ti)
	else: 
	    self.xlab = 'PGV (cm/s)'
	    self.filen = 'PGV'

        self.MetaPth = MetaPth 


    def HazardCurveCalc(self, IMLs, SiteData):
	"""
	Compute Hazard Curve and Plot 
	"""
	IMLs = IMLs.tolist() 
	MetaFile = self.MetaPth + '/Meta_HazardCurve_%s_%s_%s.py'%(IMR, self.SiteName, self.filen)

	if not os.path.exists( MetaFile ): 

	    Database = CyberShakeDatabase(self.rup_model_ids,self.Ti)
	    RupProbs = Database.ExtractRupProb() 
	    Sites = CyberShakeSites( SiteData )
	     
	    if self.IMR == 'CyberShake': 
		if self.SiteName not in Sites.keys(): 
		    print 'SiteName is not in the CyberShake Sites List'
		    raise ValueError 
		IMs = Database.ExtractCyberShakeIMs(self.SiteName) 
		PoE = []
		for iml in xrange( len(IMLs) ):
		    IML = IMLs[iml]
		    PoE0 = 1.0 
		    for sid in IMs.keys():
			skey = '%s'%sid
			Prob1 = 0.0
			for rid in IMs[skey].keys(): 
			    srkey = '%s,%s'%(sid,rid) 
			    rkey = '%s'%rid 
			    NL = len(IMs[skey][rkey])
			    NIML = len((np.array(IMs[skey][rkey])>=IML).nonzero()[0])
			    ProbIMTlgIML = NIML*1.0/NL 
			    RupProb = RupProbs[srkey] 
			    Prob1 += ProbIMTlgIML * RupProb 
			PoE0 *= (1.0 - Prob1)
		    PoE.append( 1.0 - PoE0 )
	    
	    else: 
		# Other IMR 
		if IMR == 'BA08': 
		    RrupCalc=False 
		    RxCalc = False 
		else: 
		    RrupCalc = True 
		    RxCalc = True 

		# Site info 
		SiteFlag = 0
		if self.SiteInfo['SiteName'] in Sites.keys(): 
		    print 'Given site name is in CyberShake site list'
		    SiteName = self.SiteInfo['SiteName']
		    SiteGeom = float(Sites[SiteName]['lon']), float(Sites[SiteName]['lat']),0
		    Vs30 = float(Sites[SiteName]['Vs30_Wills_2006'])
		    Z10 = float(Sites[SiteName]['Z1.0']) * 1000   # km to m
		    Z25 = float(Sites[SiteName]['Z2.5'])
		    SiteFlag = 1
		else: 
		    try: 
			SiteGeom = self.SiteInfo['lon'], self.SiteInfo['lat'],0
		    except: 
			print 'You have to specify the loation of your site'
			raise ValueError 

		    for SiteName in Sites.keys(): 
			if SiteGeom[0] == float(Sites[SiteName]['lon']) and SiteGeom[1] == float(Sites[SiteName]['lat']): 
			    print 'The site is a CyberShake site'
			    Vs30 = float(Sites[SiteName]['Vs30_Wills_2006'])
			    Z10 = float(Sites[SiteName]['Z1.0']) * 1000   # km to m
			    Z25 = float(Sites[SiteName]['Z2.5'])
			    SiteFlag = 1
			    break
		if SiteFlag == 0:
		    try: 
			Vs30 = self.SiteInfo['Vs30'] 
			Z10 = None
			Z25 = None 
		    except: 
			print 'When using other IMR but CyberShake, you should give at least the Vs30 for the site!'
			raise ValueError 

		# Source Info (fault geometry for distance calculation)
		SourceInfo = Database.ExtractSources() 
		
		mu = {}; std = {}
		for skey in SourceInfo.keys(): 
		    flon, flat, dip, ztor, zbom, rake = SourceInfo[skey]['FaultInfo'] 
		    W = (zbom-ztor)/np.sin(dip*np.pi/180.)
		    FaultTrace = []
		    for i in xrange( len(flon) ): 
			FaultTrace.append( [flon[i],flat[i],ztor] ) 
		    Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace,ztor,zbom,dip,RrupCalc=RrupCalc,RxCalc=RxCalc)
		    mu[skey] = {}; std[skey] = {}
		    for rkey in SourceInfo[skey]['Ruptures'].keys(): 
			Mw = SourceInfo[skey]['Ruptures'][rkey]

			srkey = '%s,%s'%(skey,rkey)
			median, sigmaT, tau0, sigma0 = NGA08(IMR[:2], Mw,Rjb,Vs30, self.Ti, rake=rake, dip=dip, W=W, Ztor=ztor,Rrup=Rrup,Rx=Rx, Z10=Z10,Z25=Z25) 
			mu[skey][rkey] = np.log(median[0])
			std[skey][rkey] = np.log(sigmaT[0]) 
		
		PoE = []
		for iml in xrange( len(IMLs) ):
		    IML = IMLs[iml]
		    PoE0 = 1.0 
		    for skey in mu.keys():
			Prob1 = 0.0
			for rkey in mu[skey].keys(): 
			    srkey = '%s,%s'%(skey,rkey) 
			    UpperLimit = mu[skey][rkey]+4*std[skey][rkey]
			    mu0 = mu[skey][rkey]
			    sigma0 = std[skey][rkey]
			    if np.log(IML) > UpperLimit: 
				ProbIMTlgIML = 0.0
			    else: 
				ProbIMTlgIML = GaussianCumulative(np.log(IML),UpperLimit, mu0, sigma0 )
			    RupProb = RupProbs[srkey] 
			    Prob1 += ProbIMTlgIML * RupProb 
			PoE0 *= (1.0 - Prob1)
		    PoE.append( 1.0 - PoE0 )
                
	    header = '#metafile for site %s at %s sec, IMR=%s\n'%(self.SiteName, '%.3f'%self.Ti, self.IMR)

	    # save to MetaFile 
	    meta = dict( IMLs = IMLs, 
		         PoE = PoE,)
	    save(MetaFile, meta, header=header)


	return PoE, IMLs 
    

if __name__ == '__main__':

    opt = sys.argv[1]     # CalcHC, PlotHC

    SiteInfo = {'SiteName':'STNI',}
    ERFinfo = (35,5,3,1)
    Ti = 3.0
    ERF = 'UCERF2'
    MetaPth = './metadata'
    
    if opt == 'CalcHC':
	IMR = sys.argv[2]

	P = PSHA(SiteInfo,ERFinfo, Ti, IMR=IMR, ERF=ERF,MetaPth = MetaPth  )

	IMLs = np.arange(0.01,3,0.1)   # in (g)
	SiteData = '/Users/fengw/work/Project/CyberShake_analysis/data/Sites/cs_site_types.txt'
	P.HazardCurveCalc(IMLs,SiteData)
    
    if opt == 'PlotHC':
	
	import matplotlib.pyplot as plt 
	# Plot Comparison of Hazard Curves (with CyberShake)
	IMRs = ['CyberShake','CB08','BA08','CY08','AS08']
	clrs = ['k','#FFC0CB','b','g','#FFA500']
	lss = ['-','--','--','--','--'] 
	marker=['o','','','',''] 

	P = PSHA(SiteInfo,ERFinfo, Ti, IMR='CyberShake', ERF=ERF, MetaPth = MetaPth )
	xlab = P.xlab 
	pfmt = 'png'
	fig = plt.figure(1) 
	ax = fig.add_subplot( 111 ) 
	for i in xrange( len(IMRs) ):
	    IMR = IMRs[i]
	    MetaFile = P.MetaPth + '/Meta_HazardCurve_%s_%s_%s.py'%(IMR, P.SiteName, P.filen)
	    
	    meta = load(MetaFile)
	    IMLs = meta.IMLs
	    PoE = meta.PoE 
	    ax.semilogy(IMLs,PoE,c = clrs[i], ls=lss[i], marker=marker[i],mfc='None',mec=clrs[0],label=IMRs[i])
	lg = plt.legend( title='IMRs')
	lg.draw_frame(False)
	ax.set_xlim([0,2.0])
	ax.set_ylim([1.e-6,1])
	
	plt.grid(True)
	plt.grid(b=True,which='minor')

	ax.set_xlabel( xlab ) 
	ax.set_ylabel( 'Probability Rate (1/yr)' )
	ax.set_title( 'Hazard Curve for %s'%P.SiteName )
	fig.savefig( './plots/Meta_HazardCurves_IMRs_%s_%s.%s'%(P.SiteName, P.filen, pfmt), format=pfmt )

	    

