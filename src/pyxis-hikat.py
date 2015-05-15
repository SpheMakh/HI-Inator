## HI Simulations pipeline
## github.com:SpheMakh/Ceiling-KAT_Apps [branch hikat]
## sphe sphemakh@gmail.com

# import some python essentials
import os
import sys
import math
import numpy as np
import tempfile
import scipy.ndimage as nd

# some useful constants
PI = math.pi
FWHM = math.sqrt(math.log(256))

# import pyxis essentials
import Pyxis
from Pyxis.ModSupport import *
import im
import mqt
import ms
import lsm

# Use simms to make simulatated measurement sets simms ( https://github.com/SpheMakh/simms.git )
import simms

# import non standard package
import pyfits
import Tigger
import galProps
import pyrap.measures 
dm = pyrap.measures.measures()
import pyrap.quanta as dq


import im.lwimager 

def simsky(msname='$MS',skymodel='$LSM',
           gain_err=None,pbeam=False,
           pointing_err=None,addnoise=True,sefd=551,
           noise=0,scalenoise=1,column='$COLUMN',
           freq_chunks='$FREQ_CHUNKS',
           **kw):
    """ simulate sky into MS """

    msname,skymodel,column,freq_chunks = interpolate_locals('msname skymodel column freq_chunks')
    if RECENTRE:
        recentre_fits(skymodel,ra=0.0,dec=-30.0)
    ms.set_default_spectral_info() 
    im.IMAGE_CHANNELIZE = 1

    im.lwimager.predict_vis(image=skymodel,column="MODEL_DATA",padding=1.5,wprojplanes=0)
    xo.sh('rm -fr ${MS:BASE}-${im.stokes}-channel${_nchan}.img.residual') # not sure why this image is there TODO(sphe): look into this

    if addnoise:
        noise = noise or compute_vis_noise(sefd=sefd)*scalenoise
        if scalenoise!=1:
            info('Visibility noise after scaling is %.3g mJy'%(noise*1e3))
        simnoise(noise=noise,addToCol='MODEL_DATA',column=column)
    elif column!= "MODEL_DATA":
        ms.copycol(fromcol="MODEL_DATA", tocol=column)


def azishe(fitsfile="$LSM", prefix='$MS_PREFIX', nm="$NM",
           start=1, stop=1, **kw):

    makedir(v.DESTDIR)
    
    global MS_LSM_List
    fitsfile, prefix, nm = interpolate_locals("fitsfile prefix nm")
    nm = int(nm)

    do_step = lambda step: step>=float(start) and step<=float(stop) 
    get_pw = lambda fitsname: abs(pyfits.open(fitsname)[0].header['CDELT1'])

    (ra,dec),(freq0,dfreq,nchan),naxis = fitsInfo(fitsfile)
    dfreq = abs(dfreq)
    chunks = nchan//nm
    
    mss = ['%s-%04d.MS'%(prefix,d) for d in range(nm)]

    if DEFAULT_IMAGING_SETTINGS:
        # estimate resolution using longest baseline
        tab = ms.ms()
        uvmax = max( tab.getcol("UVW")[:2].sum(0)**2 )
        res = 1.22*((2.998e8/freq)/maxuv)/2.

        # choose cellsize to be sixth of resolution
        im.cellsize = "%farcsec"%(numpy.rad2deg(res/6)*3600)
        im.npix = 2048
        im.weight = "briggs"
        im.robust = 2 # this equivalent to natural weighting

    if do_step(1):
        lsms = im.argo.splitfits(fitsfile, chunks, ctype="FREQ", prefix=prefix)
        
        with open(LSMLIST,'w') as lsm_std:
            makedir(OUTDIR)
            lsm_std.write(",".join(lsms))
            
        MS_LSM_List = ["%s,%s"%(a, b) for a, b in zip(mss, lsms) ]
        x.sh('rm -f $IMAGELIST')

        def  make_and_sim(ms_lsm_comb="$MS_LSM", options={}, **kw):
            msname,lsmname = II(ms_lsm_comb).split(',')
            v.MS = msname
            v.LSM = lsmname
            _simms(**kw)
            simsky(scalenoise=1e-8, **options)
            image()
            _add(im.RESTORED_IMAGE, IMAGELIST)

        pper('MS_LSM', lambda : make_and_sim(**kw) )
        images = get_list(IMAGELIST)
        im.argo.combine_fits(images, outname=FULL_IMAGE, ctype='FREQ', keep_old=True)


def image(msname='$MS', lsmname='$LSM', remove=None, **kw):
    """ imaging/cleaning function """
   
    global NSRC
    #NSRC = make_pure_lsm()

    im.stokes = 'I'
    im.cellseize = '3arcsec'
    im.weight = 'briggs'
    im.robust = 0
    im.npix = 2048
    msname,lsmname = interpolate_locals('msname lsmname')
    v.MS = msname
    v.LSM = lsmname

    im.IMAGE_CHANNELIZE = 1
    ms.set_default_spectral_info()


    # create clean mask to speed things up; use PURE_CAT_LSM
    #tf = tempfile.NamedTemporaryFile(suffix='.fits',dir='.')
    #tf.flush()
    #empty_image = tf.name
    #im.argo.make_empty_image(image=empty_image,channelize=0,**kw)
    #x.sh('tigger-restore ${empty_image} $PURE_CAT_LSM ${im.MASK_IMAGE} -b 20 -f')

    #mask = im.MASK_IMAGE
    #im.argo.make_threshold_mask(input='${im.MASK_IMAGE}',output=mask,threshold=0)
    #tf.close()
    
    im.IMAGER = IMAGER
    #kw.update( {"mask" if IMAGER in ["lwimager","casa"] else "fitsmask":mask} ) 
    im.make_image(dirty=False,restore=True,restored_image=RESTORED,restore_lsm=False,fitbeam=True,**kw)

    if remove:
        remove_channels_from_extremes(II(RESTORED),channels=remove)

        
def remove_channels_from_extremes(fitsname,channels=None):
    """ 
     lwimager,casa do some funny business on the first/last channel. 
     These channels will have to be removed before running the source finder;
    """
    hdu = pyfits.open(fitsname)
    data = hdu[0].data
    hdr = hdu[0].header
    # find freq axis
    axis = 0
    naxis = hdr['NAXIS']

    for i in range(1,naxis+1):
        if hdr['CTYPE%d'%i].startswith('FREQ'):
            axis = i
    chans = range(hdr['NAXIS%d'%axis])
    info('>>> Removing channels $channels from $fitsname')

    for chan in channels:
        chans.pop(chan)
    N = len(chans)
    new_fits = [slice(None)]*naxis
    new_fits[naxis-axis] = chans
    hdu[0].data = data[new_fits]

    # update hdr
    crval = hdr['CRVAL%d'%axis]
    cdelt = hdr['CDELT%d'%axis]
    crpix = hdr['CRPIX%d'%axis]
    freq0 = (crval+cdelt) - cdelt*crpix
    freq_c = freq0 + N/2*cdelt

    # use the middle channel as the reference
    hdr['CRVAL%d'%axis] = freq_c
    hdr['CRPIX%d'%axis] = N/2
    hdu[0].header = hdr
    hdu.writeto(fitsname,clobber=True)
    hdu.close()


def extract_from_skymodel(lsmname='$LSM',addnoise=True,withnoise='$WITHNOISE',checkfirst=True,**kw):
    """ run source finder on perfect cube """

    lsmname,withnoise = interpolate_locals('lsmname withnoise')
    if addnoise:
        if checkfirst:
            _min = pyfits.open(lsmname)[0].data.min()
            if _min<0:
                warn('File seems to have noise (i.e, it has negative pixels),\
                      not adding noise. Rerun with "checkfisrt=False" to force addition of noise ')
            else:
                makedir(v.DESTDIR)
                addnoise_fits(lsmname,withnoise,1e-9)

    lsm.sofia_search(withnoise,makeplot=True,**kw)

    return withnoise


def addnoise_fits(fitsname,outname,rms,correlated=False,clobber=True):
    """ add noise to fits file """

    im.argo.swap_stokes_freq(fitsname)
    hdu = pyfits.open(fitsname)
    hdr0 = hdu[0].header
    data = hdu[0].data
    noise =  np.random.normal(0,rms,data.shape) 

    if correlated:
        " do something "

    hdu[0].data = data + noise
    # Get lwimager to create a proper fits header for us
    tf = tempfile.NamedTemporaryFile(suffix='.fits',dir='.')
    tf.flush()
    tfname = tf.name
    im.argo.make_empty_image(image=tfname)
    im.IMAGER = IMAGER # Rest imager name
    hdr = pyfits.open(tfname)[0].header
    hdr['CRPIX1'] = hdr0['CRPIX1']
    hdr['CRPIX2'] = hdr0['CRPIX2']
    hdr['CDELT1'] = hdr0['CDELT1']
    hdr['CDELT2'] = hdr0['CDELT2']
    hdu[0].header = hdr
    tf.close()
    hdu.writeto(outname,clobber=clobber)
    info('Added noise to fits file, output is at $outname')
    hdu.close()

    
def _simms(msname='$MS',observatory='$OBSERVATORY',antennas='$ANTENNAS',
           synthesis="$SCAN",start_time=None,
           dtime=60,freq0='$FREQ0',dfreq='10MHz',
           nchan=2,direction='$DIRECTION',
           fromfits=True,fitsfile='$LSM',
           **kw):
    """ create simulated measurement set """

    msname,observatory,antennas,direction,freq0,fitsfile,synthesis = \
        interpolate_locals('msname observatory antennas direction freq0 fitsfile synthesis')
    if MS_REDO is False and exists(msname):
        v.MS = msname
        return

    makedir(MSOUT)

    if fromfits:
        (ra, dec),(freq0, dfreq, nchan), naxis = fitsInfo(fitsfile)

        direction="J2000,%fdeg,%fdeg"%(ra, dec)
        dfreq = abs(dfreq)
        im.argo.swap_stokes_freq(fitsfile)
        x.sh('fitstool.py -E CDELT3=$dfreq $fitsfile') 
        freq0 += dfreq/2

        if SPWIDS > 1:
            bw = nchan*dfreq
            freq0 = [f0-dfreq/2 for f0 in np.arange(freq0,freq0+bw,bw/2)]
            _nchan = nchan/SPWIDS if nchan%SPWIDS==0 else (nchan +(SPWIDS-1)) /SPWIDS
            nchan = [_nchan]*SPWIDS
        dfreq = [dfreq]*SPWIDS
                
    pos_type = "casa" if os.isdir(antennas) else "ascii"

    info('Making simulated MS: $msname ...')
    msname = simms.create_empty_ms(msname=msname, tel=observatory, direction=direction, 
                pos=antennas, pos_type=pos_type,  nbands=SPWIDS,
                synthesis=float(synthesis), dtime=dtime, freq0=freq0, 
                dfreq=dfreq, nchan=nchan, outdir=OUTDIR, **kw)

    v.MS = msname


def recentre_fits(fitsname, ra, dec):
    """ recentre FITS file """
    with pyfits.open(fitsname) as hdu:
        hdr = hdu[0].header
        hdr['CRVAL1'] = ra
        hdr['CRVAL2'] = dec
        hdu.writeto(fitsname,clobber=True)


def make_pure_lsm():
    makedir(DESTDIR)
    Model = galProps.Model(PURE_CAT,LSM)
    model = Model.load()
    model.writeto(PURE_CAT_LSM, overwrite=True)
    return model.nsrcs

def compute_vis_noise (sefd):
    """Computes nominal per-visibility noise"""

    sefd = sefd or SEFD
    tab = ms.ms()
    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")

    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID, 0]
    wavelength = 300e+6/freq0
    bw = spwtab.getcol("CHAN_WIDTH")[ms.SPWID, 0]
    dt = tab.getcol("EXPOSURE", 0, 1)[0]
    dtf = (tab.getcol("TIME", tab.nrows()-1, 1)-tab.getcol("TIME", 0, 1))[0]

    # close tables properly, else the calls below will hang waiting for a lock...
    tab.close()
    spwtab.close()

    info(">>> $MS freq %.2f MHz (lambda=%.2fm), bandwidth %.2g kHz, %.2fs integrations, %.2fh synthesis"%(freq0*1e-6, wavelength, bw*1e-3, dt, dtf/3600))
    noise = sefd/math.sqrt(2*bw*dt)
    info(">>> SEFD of %.2f Jy gives per-visibility noise of %.2f mJy"%(sefd, noise*1000))

    return noise 


def fitsInfo(fits):
    """ Returs FITS image basic info """

    hdr = pyfits.open(fits)[0].header
    ra = hdr['CRVAL1'] 
    dec = hdr['CRVAL2']
    naxis = hdr['NAXIS']

    if naxis>3: freq_ind = 3 if hdr['CTYPE3'].startswith('FREQ') else 4
    else: 
        freq_ind = 3
        if hdr['CRTYPE3'].startswith('FREQ') is False: 
            return (ra,dec), (False, False, False) , naxis

    nchan = hdr['NAXIS%d'%freq_ind]
    dfreq = hdr['CDELT%d'%freq_ind]
    freq0 = hdr['CRVAL%d'%freq_ind] + hdr['CRPIX%d'%freq_ind]*dfreq

    return (ra, dec),(freq0, dfreq, nchan), naxis


def simnoise (noise=0, rowchunk=100000, addToCol=None, scale_noise=1.0, column='MODEL_DATA'):
    """ Simulate noise into an MS """

    spwtab = ms.ms(subtable="SPECTRAL_WINDOW")
    freq0 = spwtab.getcol("CHAN_FREQ")[ms.SPWID,0]/1e6

    tab = ms.msw()
    dshape = list(tab.getcol('DATA').shape)
    nrows = dshape[0]

    noise = noise or compute_vis_noise()

    if addToCol: colData = tab.getcol(addToCol)

    for row0 in range(0,nrows,rowchunk):
        nr = min(rowchunk,nrows-row0)
        dshape[0] = nr
        data = noise*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape)) * scale_noise

        if addToCol: 
            data+=colData[row0:(row0+nr)]
            info(" $addToCol + noise --> $column (rows $row0 to %d)"%(row0+nr-1))
        else : info("Adding noise to $column (rows $row0 to %d)"%(row0+nr-1))

        tab.putcol(column,data,row0,nr)
    tab.close() 


def _add(addfile, filename='$MSLIST'):
  """ Keeps track of Files when running multiple threads.
    The files names are stored into a file which can be specified 
    via filename. The files can then be rertieved as a python list using the function get_filelist().
  """

  addfile, filename = interpolate_locals('addfile filename')

  try :
    ms_std = open(filename,"r")
    line = ms_std.readline().split("\n")[0]
    ms_std.close()
  except IOError: line = ''

  file_std = open(filename,"w")
  line+= ",%s"%addfile if len(line.split()) > 0 else addfile
  file_std.write(line)
  file_std.close()
  info("New MS list : $line")

document_globals(_add,"MSLIST")


def get_list(filename='$MSLIST'):
    """ See help for _add"""

    filename = interpolate_locals('filename')
    if not exists(filename):
        return False

    with open(filename) as ms_std:
      ms_std = open(filename)
      mslist = ms_std.readline().split('\n')[0].split(',')
    info('Found %d files.'%(len(mslist)))

    return mslist
