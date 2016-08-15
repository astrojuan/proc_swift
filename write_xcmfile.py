def write_xcmfile(xcm_array, filename="ec_srcbin10.xcm"):
    """
    for an xcm string array output from fit_ecsrcbin10, outputs the xcm string array into an xspec command file
    """
    import os
    if os.path.isfile(filename):
        print "%s exists" % filename
        ans=raw_input('Overwrite [y/n]? ')
        if ans.strip().lower()=='n':
            print "%s not overwritten; Returning" % filename
            return output, xcm
    print "Writing File %s" % filename
    f=open(filename,'w')
    for i in xcm_array:
        f.write(i+"\n")
    f.close()
    return
        
    
