#!/use/bin/env python 
"""
Intensity Measure Relations (IMRs)
"""
from utils import * 

# NGA calculation
from pynga import * 
from pynga.utils import * 

# CyberShake
class CyberShakeIMs:
    """
    CyberShake IM database 
    IMs
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


    def ExtractIMs(self,CS_Short_Name):
	
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
	Nsources = len(Source_IDs)
	itime = 0
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
		#print query, irup
		cursor.execute( query )
		IMvalues = cursor.fetchall() 
		for im in xrange( len(IMvalues) ): 
		    IMs[skey][rkey].append( IMvalues[im][0]*factor)     # unit conversion
	    
	    sys.stdout.write('\rIM extraction finished %.2f%%'%(itime*100.0/Nsources))
	    sys.stdout.flush()
	    sleep(0.1)
	    itime += 1

	    if self.Debug:
		# Debug the problem where # of PSA not equal to # Rupture Variations !
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
	sys.stdout.write('\n')

	return IMs


# get the GMPEs here input ERF site, then compute IMs for further use 
class NGA08IMs: 
    def __init__(self, IMR, Ti): 
	self.IMR = IMR
	self.Ti = Ti 

    def ExtractIMs(self, SourceInfo, SiteGeom, Vs30, Z25, Z10): 
	
	# Source info 
	flon, flat, dip, ztor, zbom, rake, Mw = SourceInfo
	W = (zbom-ztor)/np.sin(dip*np.pi/180.)
	FaultTrace = []
	for i in xrange( len(flon) ): 
	    FaultTrace.append( [flon[i],flat[i],ztor] ) 
	
	# Compute Distance Parameters
	if self.IMR == 'BA08': 
	    RrupCalc=False 
	    RxCalc = False 
	else: 
	    RrupCalc = True 
	    RxCalc = True 
	Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace,ztor,zbom,dip,RrupCalc=RrupCalc,RxCalc=RxCalc)

	# Compute NGA IMs
	median, sigmaT, tau0, sigma0 = NGA08(self.IMR[:2], Mw,Rjb,Vs30, self.Ti, rake=rake, dip=dip, W=W, Ztor=ztor,Rrup=Rrup,Rx=Rx, Z10=Z10,Z25=Z25) 
	Median = np.log(median[0])
	SigmaT = np.log(sigmaT[0]) 
	    
	return Median, SigmaT

