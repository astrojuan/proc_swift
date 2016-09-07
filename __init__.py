"""
proc_swift is a module to process and analyze swift xrt spectra
"""

from fit_srcbin10 import fit_srcbin10
from fit_srcbin10 import print_obsinfo
from fit_srcbin10 import get_obsid_flux
from reduce_xrt import reduce_xrt, mk_regions, mk_xcofiles, mk_arf, mk_spec
from reduce_xrt import create_xrtpipeline_script
from untar_swift import untar_ecswift
from write_xcmfile import write_xcmfile
from utils import create_xrtmkarf_script
from ds9 import ds9