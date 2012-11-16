#!/use/bin/env python 
""" 
Site Info 
Name, Location (lon/lat), Vs30, Vs30 flag (optional), Z2.5, Z1.0 
""" 

# CyberShake Sites
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


