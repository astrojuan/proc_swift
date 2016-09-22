def create_xrtmkarf_script(obsid, outdir, ProcDir='/Users/corcoran/Dropbox/Eta_Car/swift/quicklook',
                           TemplateDirectory='/software/github/ecswift/templates', pcmode=False):


    """
    Creates the script used to run xrtmkarf

    :param obsid: swift observation ID
    :param ProcDir: directory location of swift archived data
    :param TemplateDirectory: directory location of template script files
    :param pcmode: if True, then PC mode data, otherwise WT mode
    :return:
    """

    if pcmode:
        marftemplate = TemplateDirectory + "/mk_xrt_arfs_jamar_Temp_pc.csh"  # Location of XRT make arf file Template
    else:
        marftemplate = TemplateDirectory + "/mk_xrt_arfs_jamar_Temp.csh"  # Location of XRT make arf file Template
    marffile = "mk_xrt_arfs_jamar.csh"
    marfout = outdir + "/" + marffile  # Duplicate template
    # shutil.copy2(src5, dst5)
    # ARF = "mk_xrt_arfs_jamar.csh" #name of file to edit

    # Required: Replace all the numbers in the text, Run XRT Pipeline,
    # solution:

    tt = open(marftemplate, "r")  # opens xselect template command file
    a = tt.read()
    tt.close()
    b = a.replace("DESIREDDIRECTORY", ProcDir).replace("OBSID", obsid)
    arf = open(marfout, "w")
    arf.write(b)
    arf.close()
    print "Wrote xrtmkarf script to {0}".format(marfout)
    return marfout

def update_data(root, downloaddir, year='2016', month='08', clobber=False):
    """
    checks the Swift data/obs directory for the specified month and year, and
    downloads observation data for those observations not yet downloaded

    :param root: swift observation root (for example "00081901" for eta car data from 2016
    :param downloaddir: local directory to receive the data (ex: /Users/corcoran/research/ETA_CAR/Swift/2016)
    :param year: observation year (ex: '2016')
    :param month: observation month (2 digits, ex: '08')
    :param clobber: if True will download the data even if it currently exists
    :return: returns list of downloaded observations
    """
    import ftplib
    import os
    import shutil

    ftp = ftplib.FTP("heasarc.gsfc.nasa.gov")
    ftpdir = 'swift/data/obs/'+year.strip()+'_'+month.strip()
    ftp.login('anonymous')
    ftp.cwd(ftpdir)
    available_obs = [x for x in ftp.nlst() if root in x]
    ftp.close()

    # get list of already downloaded data
    downloaded = [x for x in os.listdir(downloaddir) if ((root in x) and ('.csh' not in x))]

    newobs=[]

    cwd = os.getcwd()

    for o in available_obs:
        if o in downloaded:
            if not clobber:
                print o + ' already exists in '+downloaddir
            else:
                newobs.append(o)
                oldobs = downloaddir+'/'+obs
                print 'Clobbering  '+oldobs
                shutil.rmtree(oldobs)
                get_swift_obs(year, month, o, downloaddir)
        else:
            print "Downloading " + o + "from Swift archive"
            # TODO: code to wget a tar file of obs data, then untar
            # note: need to retrieve obs data using
            #     wget ftp://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2016_08/00030956199.tar for example
            newobs.append(o)
            get_swift_obs(year, month, o, downloaddir)
    os.chdir(cwd)
    return newobs

def get_swift_obs(year, month, obs, DownLoadDir):
    """
    Downloads data from the swift obs archive to the download data directory
    :param year: 4 digit year (str); '2016'
    :param month: 2 digit month (str); '08'
    :param obs: Swift observation id (str); '00030956199'
    :param DownLoadDir: location to download and untar the tar file
    :return:
    """
    import os
    tarfile=obs.strip()+'.tar'
    cwd = os.getcwd()
    os.chdir(DownLoadDir)
    cmd = 'wget ftp://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/' + year.strip() + '_' + month.strip() + '/' + tarfile
    os.system(cmd)
    os.system('tar -xvf '+tarfile)
    os.system(cwd)
    os.remove(tarfile)
    os.chdir(cwd)
    return
