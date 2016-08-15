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

