#! /usr/bin/env python

"""

"""


def call_xspec():
    import xspec
    print "starting fit_srcbin10"
    import glob
    import pyfits
    import numpy as np
    from astropy.time import Time
    import xspec
    xspec.AllData.clear()
    xspec.Xset.chatter = 0
    xspec.FitManager.query="yes"
    xspec.Fit.query="yes"
    xspec.Fit.nIterations=100
    pcmode = False
    WorkDir = '/Users/corcoran/research/WR140/Swift/data/2016/work'
    rmfdir = '/caldb/data/swift/xrt/cpf/rmf'
    if pcmode:
        mode='pc'
    else:
        mode='wt'
    xspecdir=WorkDir+"/"+obsid.strip()+"/"+mode+"/xspec"
    cwd=os.getcwd()
    print "\n"
    os.chdir(xspecdir)
    src=xspec.Spectrum("srcbin10.pha")
    try:
        hdu=pyfits.open("src.arf")
    except:
        print "ARF src.arf not found; Returning"
        return
    arfhd=hdu[1].header
    try:
        respfile=arfhd['RESPFILE']
    except:
        print "RESPFILE keyword in srcbin10.arf not found; Returning"
        return
    try:
        rmffile=glob.glob(rmfdir+'/'+respfile)[0]
    except:
        print "Response file %s does not exist; Returning" % (rmfdir+'/'+respfile)
        return
    xspec.AllData.clear()
    mostring = "wabs*apec + wabs*(apec + gaussian)"
    print "test1.55"
    m = xspec.Model(mostring)
    return

if __name__ == "__main__":
    import os
    import argparse
    parser = argparse.ArgumentParser(description='see if the import of xspec causes a seg fault')
    parser.add_argument("obsid", type=str, help="Swift XRT obsid to process")
    args=parser.parse_args()
    obsid=args.obsid
    call_xspec()


