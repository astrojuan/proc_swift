import glob
import os
from astropy.io import fits as pyfits

def mk_xcofiles(obsid, mode, ProcDir='/Users/corcoran/research/WR140/Swift/data/2016'):
    """
    Creates xselect time filter  and pha extraction command files
    to be run in the work/OBSID/xsel directory for a given obsid; writes
    the .xco files in
    directory where the xselect command
    <ProcDir>/work/<OBSID>/<MODE>/xsel

    :param obsid:
    :param mode:
    :param XselectDir:
    :return:
    """

    obs = obsid.strip()

    curdir=os.getcwd()
    XselectDir = "{0}/work/{1}/{2}/xsel".format(ProcDir, obsid, mode)
    print "Changing Directory to {0}".format(XselectDir)
    os.chdir(XselectDir)

    # Find the event file
    evt = "sw{0}x{1}*po_cl.evt".format(obs,mode.strip())
    try:
        evtfile = glob.glob(evt)[0]
    except IndexError, errmsg:
        print errmsg
        print "returning"
        status = -9
        return status

    # create list of commands for time filter file
    xcofile = open('xselect{0}_timefilt.xco'.format(mode.strip()),mode='w')
    cmd = []
    cmd.append("xsel_session")
    cmd.append("$rm src.lc")
    cmd.append("$rm xtimefilt.curs_gti")
    cmd.append("read event")
    cmd.append(".")
    cmd.append("{0}".format(evtfile))
    cmd.append("yes")
    cmd.append("set device /xw")
    cmd.append("set bin 100")
    cmd.append("extract curve")
    cmd.append("save curve src")
    cmd.append("filter time cursor")
    cmd.append("re y 0 10")
    cmd.append("quit")
    cmd.append("save time cursor")
    cmd.append("xtimefilt_{0}".format(mode)) # writes file xtimefilt.curs_gti
    cmd.append("quit")
    cmd.append("no")
    for c in cmd:
        xcofile.write(c+'\n')
    xcofile.write('\n')
    xcofile.close()

    # Next create xco file to extract spectra
    xcofile = open('xselect{0}_pha.xco'.format(mode.strip()),mode='w')
    cmd=[]
    cmd.append("xsel_session")
    cmd.append("$rm src.pha")
    cmd.append("$rm bkg.pha")
    cmd.append("xsel")
    cmd.append("read event")
    cmd.append(".")
    cmd.append("{0}".format(evtfile))
    cmd.append("yes")
    cmd.append("set device /xw")
    xtfilt = "xtimefilt_{0}.curs_gti".format(mode)
    cmd.append('filter time file "{0}"'.format(xtfilt))
    cmd.append("filter region src_{0}.reg".format(mode))
    cmd.append("extract spec")
    cmd.append("save spec src")
    cmd.append("filter region bkg_{0}.reg".format(mode))
    cmd.append("extract spec")
    cmd.append("save spec bkg")
    cmd.append("quit")
    cmd.append("no")

    for c in cmd:
        xcofile.write(c+'\n')
    xcofile.write('\n')
    xcofile.close()

    # Done
    print "Changing directory to {0}".format(curdir)
    os.chdir(curdir)
    status = 0
    return status



def mk_regions(obsid, mode, ProcDir = '/Users/corcoran/research/WR140/Swift/data/2016',
               RegionDir='/Users/corcoran/research/WR140/Swift/data/2016/work/regions',
               evtfile=''):
    """
    displays the cleaned events file (of the form sw<OBSID>xpc*po_cl.evt or sw<OBSID>xwt*po_cl.evt)
    and allows the user to display and adjust the region file

    uses a template region file stored in RegionDir; saves region to
        ProcDir/work/OBSID/xsel

    Note that region file templates must be named <something>_<src|bkg>_<pc|wt>.reg for source, bkg, regions
       in pc or wt mode


    :param obsid: observation id of form '00081901001'
    :param mode: either 'pc' or 'wt' for photon counting or windowed timing
    :param ProcDir: root processing directory
    :param RegionDir: location of region files
    :return:
    """
    from ds9 import ds9
    import glob
    import os
    if not evtfile:
        try:
            fname = ProcDir + '/work/' + obsid + '/sw' + obsid + 'x' + mode.strip().lower() + '*po_cl.evt'
            evt = glob.glob(fname)[0]
        except:
            print "{0} Not Found; returning".format(fname)
            return
    else:
        evt=evtfile.strip()
    d = ds9('Reduce_XRT')
    #
    # first do the source regions
    #
    d.set('file '+evt)
    d.set('smooth')
    s = RegionDir.strip()+'/*src*'+mode.strip()+'.reg'
    try:
        sregion = glob.glob(s)[0]
    except Exception, errmsg:
        print "mk_regions: Problem getting region file like {0}".format(s)
        print errmsg
    print "Loading region {0}".format(sregion)
    d.set('regions load '+sregion)
    raw_input("\nAdjust source region if necessary; when ready, press <CR> when Done to save them ...")
    # source in frame 1
    d.set('frame 1')
    src_reg_file_dir = ProcDir+'/work/'+obsid.strip()+'/xsel/'
    if not os.path.isdir(src_reg_file_dir):
        os.makedirs(src_reg_file_dir)
    src_reg_file = src_reg_file_dir+'/src_'+mode.strip().lower()+'.reg'
    d.set('regions save '+src_reg_file)
    #
    # now background region
    #
    b = RegionDir.strip()+'/*bkg*'+mode.strip()+'.reg'
    bregion = glob.glob(b)[0]
    d.set('file '+evt)
    d.set('smooth')
    d.set('regions load '+bregion)
    raw_input("\nAdjust background region if necessary; when ready, press <CR> when Done to save them ...")
    bkg_reg_file = ProcDir+'/work/'+obsid.strip()+'/xsel/bkg_'+mode.strip().lower()+'.reg'
    d.set('regions save '+bkg_reg_file)
    d.set('exit')
    return

def mk_arf(obsid, mode,
           ProcDir='/Users/corcoran/research/WR140/Swift/data/2016'):
    obs=obsid.strip()
    xspecdir = "{0}/work/{1}/xspec".format(ProcDir, obs)
    os.chdir(xspecdir)
    exp = "{0}/work/{1}/sw{1}x{2}*po_ex.img".format(ProcDir,obs,mode)
    expofile = glob.glob(exp)[0]
    cmd = "xrtmkarf outfile={0}/src.arf " \
          "psfflag='yes' phafile = src.pha srcx = -1  srcy = -1 " \
          "clobber=yes expofile = {1}".format(xspecdir,expofile)
    arffile = open("mkarf.csh",mode="w")
    arffile.write(cmd)
    arffile.write("\n")
    arffile.close()
    os.system("source mkarf.csh")
    status = 0
    return status

def create_xrtpipeline_script(obsid, mode='wt',
                              ProcDir = '/Users/corcoran/research/WR140/Swift/data/2016',
                              clobber=clobber):
    """

    creates the xrtpipeline processing script run_xrtpipeline_<OBSID>.csh
    for the specified OBSID in ProcDir

    :param obsid: observation id (eg 00031251100)
    :param mode: either 'wt' or 'pc'
    :param ProcDir: processing directory
    :return: returns name of processing script
    """
    # get ra_obj, dec_obj from header of uf.evt file
    obs = obsid.strip()
    evtname = '{0}/{1}/xrt/event/sw{1}x{2}*po_uf.evt*'.format(ProcDir,obs, mode)
    try:
        evtfile = glob.glob(evtname)[0]
    except Exception, errmsg:
        print errmsg
        print evtname
        print "returning"
        return
    hdu = pyfits.open(evtfile)
    ra_obj = hdu[0].header['RA_OBJ']
    dec_obj = hdu[0].header['DEC_OBJ']

    cmd = []
    if clobber:
        clob='yes'
    else:
        clob='no'
    scriptname = ProcDir+'/run_xrtpipeline_'+obsid.strip()+'.csh'
    f = open(scriptname,mode = 'w')
    cmd.append('xrtpipeline srcra="{0:.4f}" srcdec="{1:.4f}" indir={2} '
               'outdir=work/{2}/xrtfiles steminputs="{2}", clobber="{3}"'.format(ra_obj, dec_obj, obs, clob))
    cmd.append('mkdir -p work/{0}/xrtfiles'.format(obs,mode))
    cmd.append('cd work/{0}'.format(obs))
    cmd.append('mv *.* xrtfiles/. # moving xrtpipeline output to xrtfiles')
    cmd.append('mkdir -p {0}/xsel'.format(mode))
    cmd.append('mkdir -p {0}/xspec'.format(mode))
    cmd.append('\n')
    cmd.append('cd {0}/xsel'.format(mode))
    evtf = 'sw{0}x{1}*po_cl.evt'.format(obs, mode)
    cmd.append('echo "linking cleaned event file {0}"')
    cmd.append("ln -nfs `ls {0}/work/{1}/xrtfiles/{2}` .".format(ProcDir,obs,evtf))
    cmd.append('cd ../xspec')
    for c in cmd:
        f.write(c+'\n')
    f.write('\n')
    f.close()
    return scriptname

def reduce_xrt(obsid,
               ProcDir = '/Users/corcoran/research/WR140/Swift/data/2016',
               pcmode=False,
               clobber = 'True'):
    """
    For a given obsid, this routine runs 
      1) DATA REDUCTION: runs xrt pipeline on the downloaded data (in the ProcDir directory) and outputs
      to the work directory (work/OBSID/xrtfiles)
      2) DATA EXTRACTION: creates xselect command files in work/OBSID/xsel/{mode} (where {mode} is 'wt',
      'pc' thenruns xselect to time filter, region filter and extract the source and background spectra,
      and copy them to the xsel/{mode} work  directory;
      Directory structure:
      ProcDir is where the data directory from the archive (or untarred from the quicklook data) is located
      xrtpipeline outputs to ProcDir/work/obsid
      xsel is run in ProcDir/work/obsid/xsel
      xspec is run in ProcDir/work/obsid/xspec

      changes:
      20150818 MFC added flag for analysis of pcmode data (and created associated template file)
     """
    import os
    import pyfits
    import urllib

    obs = obsid.strip()
    WorkDir = ProcDir+"/work"
    print "Writing Output to %s" % WorkDir
    #
    # need to have CALDB defined for xrtpipeline to work
    #
    try:
        caldb=os.environ['CALDB']
    except:
        print "CALDB not defined; Returning"
        return # exits the script
    # print "$CALDB = %s" % caldb
    if pcmode:
        mode = 'pc'
    else:
        mode = 'wt'
    #
    # Run XRT Pipeline
    #
    XRT = create_xrtpipeline_script(obs, ProcDir=ProcDir, mode=mode, clobber = clobber)
    os.chdir(ProcDir)
    os.system("source "+ XRT)
    #
    # we now have the cleaned events files in the xrtfiles subdirectory
    # now do the data filtering using xselect
    #
    return

def mk_spec(obsid,
              ProcDir = '/Users/corcoran/research/WR140/Swift/data/2016',
              mode = 'wt'):
    """
    Creates the region and time filter files for xselect, then runs xselect to
    create the source and background spectra for the given obsid and mode
    :param obsid:
    :param ProcDir:
    :param mode:
    :return:
    """
    #
    # display source and background regions and adjust if necessary
    #
    mk_regions(obsid,mode,ProcDir=ProcDir)
    #
    # Run Xselect
    #
    obs = obsid.strip()
    mode = mode.strip()
    XselectDir = "{0}/work/{1}/{2}/xsel/".format(ProcDir, obs, mode)
    if os.path.exists(XselectDir):
        # print "Creating xselect directory {0}".format(XselectDir)
        # os.makedirs(XselectDir)
        os.chdir(XselectDir)
        status = mk_xcofiles(obsid, mode)
        if status <> 0:
            print "Error in mk_xcofiles; returning"
            return
        cmd = "xselect @xselect{0}_timefilt.xco".format(mode.strip())
        print cmd
        os.system(cmd)
        #Create time filtered, region filtered pha src & bkg files
        cmd = "xselect @xselect{0}_pha.xco".format(mode.strip())
        print cmd
        os.system(cmd)
        #
        # change to xspec directory & copy pha files
        #
        XspectDir = ProcDir+"/work/"+obsid+"/xspec/" #specifies the location of the Xspec directory
        os.chdir(XspectDir)# Changed directory to the xspec directory
        os.system("cp ../xsel/*pha .")
        print "Current working dir :",os.getcwd() #Displays current directory
        # group the pha files
        os.system("grppha src.pha srcbin10.pha comm = 'group min 10' temp = exit") #runs the grppha command file.
        #### Make Ancillary response Files

        #marf_script = create_xrtmkarf_script(obsid, XspectDir, ProcDir=ProcDir, pcmode=pcmode)
        #os.system("source "+ marf_script) #make Ancillary Response File.
        mk_arf(obsid,mode, ProcDir=ProcDir)

        #### Locate the Response file
        hdulist = pyfits.open('src.arf') #open the Fits file.
        Response = hdulist[1].header['RESPFILE']
        hdulist.close()

        print "Getting rmf file %s from remote CALDB" % Response
        swiftrmfdir="http://heasarc.gsfc.nasa.gov/FTP/caldb/data/swift/xrt/cpf/rmf"
        urllib.urlretrieve(swiftrmfdir+"/"+Response,XspectDir+"/"+Response)
    else:
        print "{0} Does not exist; returning".format(XselectDir)
    return


if __name__ == "__main__":
    obsid = '00031251110'
    mode = 'wt'
    procdir = '/Users/corcoran/research/WR140/Swift/data/2016/ql'
    reduce_xrt(obsid,
               ProcDir=procdir,
               pcmode=False)


