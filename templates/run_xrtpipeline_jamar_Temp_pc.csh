
echo 
echo NOW PROCESSING 00000000000 NUMBER        1 OUT OF       1
echo
xrtpipeline srcra="161.264775" srcdec="-59.684431" indir=00000000000 outdir=work/00000000000 steminputs="00000000000"

mkdir work/00000000000/xsel work/00000000000/xspec

cd work/00000000000/xsel
ln -nfs `ls ../sw*po_cl.evt` .
cp /Users/corcoran/Dropbox/Eta_Car/swift/quicklook/templates/*reg .

cd ../xspec
cp /Users/corcoran/Dropbox/Eta_Car/swift/quicklook/templates/ec*xcm .

cd ../../..
#
##
#


