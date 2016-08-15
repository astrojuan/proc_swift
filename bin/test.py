for i in range(len(mseq)):
    ms = mseq[i]
    try:
        ind=where(ttseq == mseq[i])[0][0]
        mjd=ttmjd[ind]
        try:
            ms = int(ms)
            ind = where(js == ms)
            mf=mflux[i]
            mfe=mfluxerr[i]
            try:
                ind=ind[0][0] # get the index from the tuple
                jf=jflux[ind]
                jfe=jfluxerr[ind]
                flux=average([mf,jf])
                diff=abs(jf-mf)
                fluxerr=sqrt(jfe**2+mfe**2+diff**2)
                tf.add_row([ms,mjd,mf,mfe,jf,jfe,mf-jf,flux,fluxerr])
                #print("{0} {1:20.4f} {2:.3e} {3:.3e} {4:.3e} {5:+.3e} {6:.3e} {7:.3e} {8:.3e}".format(ms, jdmid, mf, mfe, jf, jfe, mf-jf, flux, fluxerr))
            except:
                print "Jamar Sequence {0} Not found".format(ms)
        except:
            print "Problem converting {0} to integer".format(ms)
    except:
        print "Problem with sequence {0}".format(mseq[i])
tf.remove_row(0)
print tf
