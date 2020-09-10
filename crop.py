from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import argparse
#https://github.com/revoltek/scripts/blob/master/lib_fits.py
from lib_fits import flatten

parser = argparse.ArgumentParser(description='Crop FITS file')
parser.add_argument('image',nargs=1,help='Image to crop')
parser.add_argument('--region',dest='region',nargs=4,type=float,default=[137.307,51.546,0.6,0.4],help='Crop region (<center ra> <center dec> <width ra> <width dec>) (all in deg)')
parser.add_argument('--output',dest='output',default='crop.fits',help='Name of output file (default crop.fits)')
args = parser.parse_args()

filename = args.image[0]
ra = args.region[0]
dec = args.region[1]
ra_width = args.region[2]
dec_width = args.region[3]

#open file, extract wcs data
hdr, data = flatten(filename)
w = WCS(hdr)

#find corners
ara = ra + ra_width / 2
adec = dec + dec_width / 2
bra = ra + ra_width / 2
bdec = dec - dec_width / 2
cra = ra - ra_width / 2
cdec = dec + dec_width / 2
dra = ra - ra_width / 2
ddec = dec - dec_width / 2

#convert corners into pixel tupels
apxx = w.wcs_world2pix(ara, adec, 1)[0]
apxy = w.wcs_world2pix(ara, adec, 1)[1]
bpxx = w.wcs_world2pix(bra, bdec, 1)[0]
bpxy = w.wcs_world2pix(bra, bdec, 1)[1]
cpxx = w.wcs_world2pix(cra, cdec, 1)[0]
cpxy = w.wcs_world2pix(cra, cdec, 1)[1]
dpxx = w.wcs_world2pix(dra, ddec, 1)[0]
dpxy = w.wcs_world2pix(dra, ddec, 1)[1]

xmax = np.rint(max([apxx, bpxx, cpxx, dpxx])).astype(np.int32)
xmin = np.rint(min([apxx, bpxx, cpxx, dpxx])).astype(np.int32)
ymax = np.rint(max([apxy, bpxy, cpxy, dpxy])).astype(np.int32)
ymin = np.rint(min([apxy, bpxy, cpxy, dpxy])).astype(np.int32)

#crop image
cropdata = data[ymin:ymax,xmin:xmax]

#find new center pixel
centerpx_y = np.rint(cropdata.shape[0] / 2).astype(np.int32)
centerpx_x = np.rint(cropdata.shape[1] / 2).astype(np.int32)

#update header data
hdr['NAXIS1'] = cropdata.shape[0]
hdr['NAXIS2'] = cropdata.shape[1]
hdr['CRPIX1'] = centerpx_x
hdr['CRPIX2'] = centerpx_y
hdr['CRVAL1'] = ra
hdr['CRVAL2'] = dec

#create new cropped fits file
fits.writeto(args.output,cropdata,hdr,overwrite=True)