

import os
import xspec
import sys

def fit_srcbin10(obsid, WorkDir="/Users/corcoran/Dropbox/Eta_Car/swift/quicklook/work",
                   emin=2.0, emax=10.0, rmfdir='/caldb/data/swift/xrt/cpf/rmf', chatter=0, pcmode = False):
    print "starting fit_srcbin10"
    import glob
    import pyfits
    import numpy as np
    from astropy.time import Time
    xspec.AllData.clear()
    xspec.Xset.chatter = chatter
    xspec.FitManager.query="yes"
    xspec.Fit.query="yes"
    xspec.Fit.nIterations=100
    if pcmode:
        mode='pc'
    else:
        mode='wt'
    xspecdir=WorkDir+"/"+obsid.strip()+"/"+mode+"/xspec"
    cwd=os.getcwd()
    print "\n"
    os.chdir(xspecdir)
    # change to use unbinned data and cstat to avoid problems with low count rate binned data

    xspec.Fit.statMethod = "cstat"
    src = xspec.Spectrum("src.pha")
    xspec.Plot.setRebin(minSig=3.0,  maxBins=10.0)
    #src=xspec.Spectrum("srcbin10.pha")
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
    src.response=rmffile
    src.response.arf="src.arf"
    src.background="bkg.pha"
    src.ignore("0.0-1.0 7.9-**")
    hdulist=pyfits.open("srcbin10.pha")

    prihdr = hdulist[0].header
    """
    dateobs=prihdr['DATE-OBS']
    dateend=prihdr['DATE-END']
    t=Time(dateobs,scale='utc', format='isot')
    te=Time(dateend,scale='utc',format='isot')
    """
    mjds = prihdr['MJD-OBS']
    mjde = prihdr['MJD-OBS']+(prihdr['TSTOP']-prihdr['TSTART'])/86400.0
    t = Time(mjds,scale='utc', format='mjd')
    dateobs = t.iso.replace(' ', 'T')
    te = Time(mjde,scale='utc', format='mjd')
    tmid=(te.jd-t.jd)/2.0+t.jd
    xspec.AllData.ignore("bad")


    """
    first fit without the cflux component
    """
    mostring="wabs*apec + wabs*(apec + gaussian)"
    m=xspec.Model(mostring)
    m.wabs.nH=[1.1,0.1,0,0,5,5]
    m.apec.kT=[1.0,0.1,0.1,0.1,2,2]
    m.apec.norm=[0.1,0.1,0.0,0.0,2,2]
    m.wabs_3.nH=[10,0.1,0,0,100,100]
    m.apec_4.kT=[4.5,-0.1,2.0,2.0,6.0,6.0]
    m.apec_4.norm=[0.1,0.1,0.0,0.0,2,2]
    m.apec_4.Abundanc=0.4
    m.gaussian.LineE=[6.4,-0.1]
    m.gaussian.Sigma=[0.01,-0.01]
    m.gaussian.norm=[1e-4,1e-5,0,0,0.1,0.1]
    m.show()

    try:
        xspec.Fit.perform()
    except Exception, errmsg:
        sys.exit("{0}; Can't proceed; exiting".format(errmsg.value))
    m.show()
    xspec.Plot.xAxis = "KeV"
    xspec.Plot.device = "/xw"
    xspec.Plot.add = True

    xspec.Plot.addCommand("re y 1e-5 20")
    xspec.Plot("ldata")
    xspec.AllModels.calcFlux(str(emin)+" "+str(emax))
    oflux=src.flux[0]

    """
    Create XCM output
    """
    xcm=['data '+src.fileName]
    xcm.append('back '+src.background.fileName)
    xcm.append('resp '+src.response.rmf)
    xcm.append('arf '+src.response.arf)
    xcm.append('ignore '+src.ignoredString())
    xcm.append('model '+ m.expression)
    for i in np.arange(m.nParameters)+1:
        p=m(i).values
        parvals=str(p[0])
        for k in np.arange(len(p)-1)+1:
            parvals=parvals+', '+str(p[k])
        xcm.append(parvals)
    xcm.append('statistic '+xspec.Fit.statMethod)
    nh=m.wabs.nH.values[0]
    kt=m.apec.kT.values[0]
    norm=m.apec.norm.values[0]
    nh3=m.wabs_3.nH.values[0]
    xspec.Fit.error("2.706 6")  # param 6 is nh3
    p6 = xspec.AllModels(1)(6)
    nhmin=p6.error[0]
    nhmax=p6.error[1]
    kt4=m.apec_4.kT.values[0]
    norm4=m.apec_4.norm.values[0]
    gnorm=m.gaussian.norm.values[0]
    """
    now add cflux component, set the values to the previous best fit,
    freeze apec norms, and refit to get the cflux and its errors
    """
    mostring2="wabs*apec + cflux*wabs*(apec + gaussian)"
    m2=xspec.Model(mostring2)
    m2.wabs.nH=[nh,-0.1]
    m2.apec.kT=[kt,-0.1]
    m2.apec.norm=[norm,-0.1]
    m2.wabs_4.nH=[nh3,-0.1]
    m2.apec_5.kT=[kt4,-0.1]
    m2.apec_5.norm=[norm4, -0.0001]
    m2.apec_5.Abundanc=0.4
    m2.gaussian.LineE=[6.4,-0.1]
    m2.gaussian.Sigma=[0.01,-0.01]
    m2.gaussian.norm=[gnorm,-0.0001]
    m2.cflux.Emin=emin
    m2.cflux.Emax=emax
    try:
        xspec.Fit.perform()
        m2.show()
        cf = m2.cflux.lg10Flux.values[0]
        print "%5.1f - %5.1f keV flux from cflux = %10.4e" % (emin, emax, 10 ** cf)
        xspec.Fit.error("2.706 8")  # param 8 is cflux
        p8 = xspec.AllModels(1)(8)
        cfmin = p8.error[0]
        cfmax = p8.error[1]
        fluxplus = 10 ** cfmax - 10 ** cf
        fluxmin = 10 ** cf - 10 ** cfmin
    except Exception, errmsg:
        print str(errmsg)
        fluxplus = 0.0
        fluxmin = 0.0
    duration=te.jd-t.jd
    nhpluserr = nhmax-nh3
    nhminerr  = nh3-nhmin
    output={'obsid':obsid, 't_start':t.jd, 't_end':te.jd, 't_mid':tmid, 'duration':duration,
            'exposure':src.exposure, 'flux':oflux, 'fluxplus':fluxplus, 'fluxmin':fluxmin,
            'dateobs':dateobs, 'nh':nh3, 'nhplus':nhpluserr, 'nhmin':nhminerr, 'emin':emin, 'emax':emax}

    """
    Write xcm file; warn user if it exists and give option to overwrite
    """
    filename=xspecdir+'/srcbin10.xcm'
    if os.path.isfile(filename):
        print "File %s exists" % filename
        ans=raw_input('Overwrite [y]/n? ')
        if ans.strip().lower()=='n':
            print "File %s not overwritten; Returning" % filename
            return output, xcm
    print "Writing file %s" % filename
    f=open(filename,'w')
    for i in xcm:
        f.write(i+"\n")
    f.close()
    #print "Changing directory back to %s" % cwd
    print "\n"
    obsinfo_log = xspecdir+"/obsinfo.log"
    print_obsinfo(output, outfile=obsinfo_log)
    return output, xcm

def print_obsinfo(obsinfo,outfile='',append=False):
    """
    prints the obsinfo dictionary returned by fit_ecsrcbin10 in a formatted way on the screen
    or optionally into a file.  If file output is requested, the data can be appended onto the file
    @param obsinfo: the obsinfo dictionary returned by fit_ecsrcbin10
    @param outfile: optional;  name of file to write output to
    @param append: optional; if True, the output is appended to the named file
    @return:
    """
    output = '{obsid:11s} {t_start:11.4f} {t_end:11.4f} {t_mid:11.4f} {duration:6.4f} {exposure:5.2f} {flux:6.3e} {fluxplus:+6.3e}  ' \
             '{fluxmin:+6.3e} {dateobs:23s} {nh:5.2f} {nhplus:+5.2f} {nhmin:+5.2f} {emin:4.2f} {emax:4.2f}'.format(**obsinfo)
    print output
    if outfile:
        if append:
            print "Appending to file {}".format(outfile)
            file=open(outfile,'a')
            file.write(output)
            file.write("\n")
            file.close()
        else:
            overwrite=True
            if os.path.isfile(outfile):
                print "{} exists".format(outfile)
                ans = raw_input('Overwrite [y]/n? ')
                if ans.strip().lower() == 'n':
                    overwrite = False
            if overwrite:
                print "Writing to file {}".format(outfile)
                file = open(outfile, 'w')
                file.write(output)
                file.write("\n")
                file.close()
    return


def get_obsid_flux(obsid, emin=2.0, emax=10.0,
                   workdir="/Users/corcoran/research/ETA_CAR/Swift/2014/work",
                   outfile='', append=False):
    """
    returns the flux for a given eta car obsid

    @param obsid:
    @param emin:
    @param emax:
    @param workdir:
    @param outfile:
    @param append:
    @return:
    """
    obsinfo, xcm = fit_ecsrcbin10(obsid, WorkDir=workdir, emin=emin, emax=emax, chatter=0)
    if not outfile:
        outfile = workdir + "/" + obsid.strip() + "/xspec/obsinfo.log"
    print_obsinfo(obsinfo, outfile=outfile, append=append)
    return

def run_fit(obsid, mode, WorkDir='/Users/corcoran/research/WR140/Swift/data/2016/work',
            emin=2.0, emax=10.0, rmfdir='/caldb/data/swift/xrt/cpf/rmf'):
    """
    for a single obsid or a list of obsids fit the swift spectra using fit_srcbin10 and
    save the output
    :param obsid:
    :param mode:
    :param WorkDir:
    :param emin:
    :param emax:
    :param rmfdir:
    :return:
    """
    import numpy as np
    if np.isscalar(obsid):
        obsid = [obsid] # create an array if a scalar value passed in
    for obs in obsid:
        obsinfo, xcm = fit_srcbin10(obsid, WorkDir=WorkDir, emin=emin, emax=emax, rmfdir=rmfdir)
        obsinfo_log=WorkDir+"/"+obsid.strip()+"/xspec/obsinfo.log"
        print_obsinfo(obsinfo, outfile=obsinfo_log)
        resp = raw_input("Press <CR> to Finish :")
    return

def check_response(question):
    resp = raw_input(question+" [y]/n :") or "y"
    if resp.strip().lower()[0] == "y":
        response=True
    else:
        response=False
    return response

if __name__ == "__main__":
    #from proc_swift import fit_srcbin10, print_obsinfo
    # import argparse
    # parser = argparse.ArgumentParser(description='Fit XRT spectrum for swift data;  '
    #                                              'the analysis results are in the ProcDir/work/OBSID/xspec subdirectory')
    # parser.add_argument("obsid", type=str, help="Swift XRT obsid to process")
    # parser.add_argument("--rmfdir", help="rmf directory (default /caldb/data/swift/xrt/cpf/rmf)", type=str, default="/caldb/data/swift/xrt/cpf/rmf")
    # parser.add_argument("--ProcDir", help="Directory where Processed data is located (default = /Users/corcoran/Dropbox/Eta_Car/swift/quicklook)", type=str, default="/Users/corcoran/Dropbox/Eta_Car/swift/quicklook")
    # args=parser.parse_args()
    # obsid=args.obsid
    # ProcDir=args.ProcDir
    # rmfdir=args.rmfdir
    # WorkDir = ProcDir+'/work'
    # #fit spectrum
    # obsinfo, xcm = fit_srcbin10(obsid, WorkDir=WorkDir, emin=2.0, emax=10.0, rmfdir=rmfdir)
    # # record info
    # obsinfo_log=WorkDir+"/"+obsid.strip()+"/xspec/obsinfo.log"
    # print_obsinfo(obsinfo, outfile=obsinfo_log)
    # resp = raw_input("Press <CR> to Finish :")
    obsid = '00031251003'
    WorkDir = '/Users/corcoran/research/WR140/Swift/2009/work'
    fit_srcbin10(obsid, WorkDir=WorkDir,
                 emin=2.0, emax=10.0, rmfdir='/caldb/data/swift/xrt/cpf/rmf', chatter=0, pcmode=True)
    print "starting fit_srcbin10"
