def untar_ecswift(obsid, tardir="/Users/corcoran/Downloads", ProcDir="/Users/corcoran/Dropbox/Eta_Car/swift/quicklook"):
    """
    this moves the tar from the tardir to ProcDir and untars it
    """
    import tarfile
    import os, glob
    import shutil
    from pprint import pprint
    
    #Required: Locate the tar file and copy the file into the desired directory.
    #Solution:
    homedir=os.getcwd()
    #os.chdir(tardir)# Changed directory to the  directory containing the tar file
    tarname="sw"+obsid+"*.tar*" # locates files with names like sw00091911037.003.tar.gz
    try :
        TarFilename = glob.glob(tardir+"/"+tarname)[-1] 
    except: 
        print "Cannot find %s; Returning" % (tarname)
        #print "The current working directory contains the following files: "
        #pprint(os.listdir(tardir))
        return
    Tname=TarFilename[TarFilename.rfind('/')+1:] # strip off directory path
    #print "Current directory : %s" % os.getcwd()
    print "Moving %s to %s " % (TarFilename, ProcDir)
    try:
        shutil.move(TarFilename,ProcDir+"/.") #moves the tar file to ProcDir
    except:
        print "Could not move %s from %s to %s/." % (TarFilename, tardir, ProcDir)
        return
    # untar the file 
    os.chdir(ProcDir)
    theTarFile = Tname # tar file to extract
    extractTarPath = '.' # tar file path to extract
    try:
        tarfile.is_tarfile(theTarFile)
    except:
        print theTarFile + " is not a tarfile; returning"
        return
    tfile = tarfile.open(theTarFile) # open the tar file
    # list all contents
    print "tar file contents:"
    print tfile.list(verbose=False)
    # extract all contents
    tfile.extractall(extractTarPath)
    print "Untar complete!"
    os.chdir(homedir)
    return
