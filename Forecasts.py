#!/usr/bin/env python 
""" 
Class of Earthquake Rupture Forecasting database 
""" 
from utils import *


# Model 1: UCERF2 (for Southern California)
class UCERF2Database: 
    
    def __init__(self, erf_id=35, rup_var_id=3):
	hostname = 'focal.usc.edu'  # host                            
	username = 'cybershk_ro'    # usrname                         
	password = 'CyberShake2007' # password                        
	database = 'CyberShake'     # Database name 
	try:
	    self.db = mdb.connect( host = hostname, user = username, \
		    passwd = password, db = database )
	except:
	    self.db = None
        
	self.erf_id = erf_id 
	self.rup_var_id = rup_var_id


    def ExtractRupProb(self): 
	cursor = self.db.cursor()
	query = 'Select Source_ID, Rupture_ID, Prob from Ruptures where ERF_ID=%s'%self.erf_id 
	cursor.execute( query ) 
	rows = cursor.fetchall()
	RupProb = {}
	for irow in xrange( len( rows ) ):
	    sid = rows[irow][0]
	    rid = rows[irow][1]
	    srkey = '%s,%s'%(sid,rid)
	    prob = rows[irow][2]
	    RupProb[srkey] = prob 
        cursor.close() 
	return RupProb


    def ExtractSource0(self): 
	cursor = self.db.cursor()
	query = 'Select Source_ID from Ruptures where ERF_ID=%s and Rupture_ID=0'%(self.erf_id)
	cursor.execute( query ) 
	Sids = cursor.fetchall() 
	SourceInfo0 = {}
	for isid in xrange( len(Sids) ): 
	    skey = '%s'%isid 
	    query = 'Select Rupture_ID from Ruptures where ERF_ID=%s and Source_ID=%s'%(self.erf_id,skey)
	    cursor.execute( query ) 
	    rids = cursor.fetchall() 
	    SourceInfo0[skey] = len(rids)
	
	return SourceInfo0 




    def ExtractRuptures(self,skey,rkey): 
	# one Source's rupture information (detailed)
	
	cursor = self.db.cursor()
	query = 'Select Mag from Ruptures where ERF_ID=%s and Source_ID=%s and Rupture_ID=%s'%(self.erf_id,skey,rkey)
	cursor.execute( query ) 
	Mw = cursor.fetchone()[0]

	# Fault Trace (find another way to get the fautl trace based on deformation model, not all discretized points)
	query = "select * from %s where ERF_ID = %s and Source_ID = %s and Rupture_ID=%s"%('Points',self.erf_id,skey,rkey)
	cursor.execute( query )       # run query
	row_point = np.array(cursor.fetchall())
	Ztor = row_point[0,6]
	Zbom = row_point[-1,6] 
	rake = np.mean( row_point[:,7] )
	dip = np.mean( row_point[:,8] )
	index_surf = (row_point[:,6] == Ztor).nonzero()[0]
	lon = row_point[index_surf,5].tolist()
	lat = row_point[index_surf,4].tolist()
	SourceInfo = [lon,lat,dip,Ztor,Zbom,rake,Mw] 
	cursor.close()
        
        return SourceInfo 


    def ExtractRuptureVariation(self):
	"""
	Extract Rupture Variations (just for this UCERF2 database) 
	""" 
	# hypocenter and slip distribution informations
	cursor = self.db.cursor()
	cursor.close()



