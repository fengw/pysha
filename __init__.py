#!/usr/bin/env python
"""
Test of the package
""" 

import Forecasts
import IMRs
import SHA 
import Sites
from utils import * 

wrk = '/Users/fengw/local/pylib/pysha'
datapth = os.path.join(wrk,'data/')
MetaPth = os.path.join(wrk,'metadata')
plotpth = os.path.join( wrk, 'plots' ) 

SiteInfo = {'SiteName':'STNI',}
SiteInfo = {'SiteName':'s758',}
SiteInfo = {'SiteName':'s115',}

ERF = 'UCERF2'
SourceType = 'NonPoisson'
ERFinfo = (35,5,3,1)   # will change depend on different IMRs and different ERFs
Ti = 3.0
IMTdict = {'-1.0':'PGA', '-2.0':'PGV', 'T':'PSA'}

if Ti != -2.0: 
    if Ti == -1.0:
	filen = 'PGA'
	xlab = 'PGA (g)'
    else: 
	filen = 'PSA%s'%('%.3f'%Ti)
	xlab = '%s sec PSA (g)'%('%.3f'%Ti)
else: 
    xlab = 'PGV (cm/s)'
    self.filen = 'PGV'

pfmt = 'png'

if __name__ == '__main__':

    import os, sys 
    opt = sys.argv[1]     # CalcHC, PlotHC, Disaggregation

    
    if opt == 'CalcHC':
	#IMRs = ['CyberShake','CB08','BA08','CY08','AS08']
	IMR = sys.argv[2]

	IMLs = np.arange(0.01,3,0.2)   # in (g)
	P = SHA.PSHA(ERF=ERF, IMR=IMR, SourceType=SourceType, ERFinfo=ERFinfo, Ti=Ti)

	SiteData = datapth + 'Sites/cs_site_types.txt'
	P.HazardCurveCalc(IMLs,SiteInfo,SiteData,MetaPth = MetaPth)
    

    if opt == 'PlotHC':

	IMR = sys.argv[2] 
	SiteName = SiteInfo['SiteName'] 
	plotpth = plotpth + '/%s'%SiteName 
	if not os.path.exists(plotpth): 
	    os.mkdir( plotpth )

	# Plot Comparison of Hazard Curves (with CyberShake)
	IMR_all = ['CyberShake','CB08','BA08','CY08','AS08']
	clrs = ['k','#FFC0CB','b','g','#FFA500']
	lss = ['-','--','--','--','--'] 
	marker=['o','','','',''] 
	
	if IMR == 'IMRs': 
	    IMRmodel = IMR_all
	    itmp = 0
	else: 
	    IMRmodel = [IMR,]
	    for i in xrange( len(IMR_all) ): 
		if IMR == IMR_all[i]: 
		    itmp = i
		    break
	    
	import matplotlib.pyplot as plt 
	fig = plt.figure(1) 
	ax = fig.add_subplot( 111 ) 

	for i in xrange( len(IMRmodel) ):
	    IMR0 = IMRmodel[i]
	    MetaFile = MetaPth + '/%s/Meta_HazardCurve_%s_%s_%s.py'%(SiteName, IMR0, SiteName, filen)
	    
	    meta = load(MetaFile)
	    IMLs = meta.IMLs
	    PoE = meta.PoE 

	    ax.loglog(IMLs,PoE,c = clrs[i+itmp], ls=lss[i+itmp], marker=marker[i+itmp],mfc='None',mec=clrs[0],label=IMR_all[i+itmp],lw=2)
	lg = plt.legend( title='IMRs' )
	lg.draw_frame(False)
	ax.set_xlim([0,2.0])
	ax.set_ylim([1.e-6,1])
	
	plt.grid(True)
	plt.grid(b=True,which='minor')

	ax.set_xlabel( xlab ) 
	ax.set_ylabel( 'Probability Rate (1/yr)' )    # this depends on the rupture probability used in UCERF2 (whether it is 1 year span or T year span)
	ax.set_title( 'Hazard Curve for %s'%SiteName )
	fig.savefig( plotpth + '/HazardCurves_%s_%s_%s.%s'%(IMR, SiteName, filen, pfmt), format=pfmt )

	    

