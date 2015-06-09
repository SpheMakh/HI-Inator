## HI Simulations pipeline
## github.com:SpheMakh/Ceiling-KAT_Apps [branch hikat]
## sphe sphemakh@gmail.com

# import some python essentials
import os
import sys
import math
import numpy as np
import tempfile
import json
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
from Simms import simms

# import non standard package
import pyfits
import Tigger
import galProps
import pyrap.measures 
dm = pyrap.measures.measures()
import pyrap.quanta as dq


import im.lwimager 


def simsky(msname='$MS', skymodel='$LSM',
           gain_err=None, pbeam=False,
           pointing_err=None, addnoise=True, sefd=551,
           noise=0, scalenoise=1, column='$COLUMN',
           freq_chunks='$FREQ_CHUNKS',
           **kw):
    """ simulate sky into MS """

    msname, skymodel, column, freq_chunks = interpolate_locals('msname skymodel column freq_chunks')

    if RECENTRE:
        ra, dec = map(np.rad2deg, [ dm.direction(DIRECTION)[m]["value"] for m in "m0","m1"])
        recentre_fits(skymodel, ra=ra, dec=dec)

    ms.set_default_spectral_info()
    im.IMAGE_CHANNELIZE = 1

    im.lwimager.predict_vis(image=skymodel, column="MODEL_DATA", padding=1.5, wprojplanes=128)
    # Clean up after lwimager
    xo.sh('rm -fr ${MS:BASE}*-channel*.img.residual') 

    if addnoise:
        noise = noise or compute_vis_noise(sefd=sefd)*scalenoise
        if scalenoise!=1:
            info('Visibility noise after scaling is %.3g mJy'%(noise*1e3))
        simnoise(noise=noise, addToCol='MODEL_DATA', column=column)

    elif column!= "MODEL_DATA":
        ms.copycol(fromcol="MODEL_DATA", tocol=column)


_ANTENNAS = {
    "meerkat": "MeerKAT64_ANTENNAS",
    "kat-7": "KAT7_ANTENNAS",
    "jvlaa": "vlaa.itrf.txt",
    "jvlab": "vlab.itrf.txt",
    "jvlac": "vlac.itrf.txt",
    "jvlad": "vlad.itrf.txt",
    "wsrt": "WSRT_ANTENNAS",
    "ska1mid254": "SKA1REF2_ANTENNAS",
    "ska1mid197": "SKA197_ANTENNAS",
}

def azishe(fitsfile="$LSM", prefix='$MS_PREFIX', nm="$NM",
            clean=True, dirty=False, psf=False,
           start=1, stop=1, config=None, addnoise=True):

    makedir(v.DESTDIR)
    
    global MS_LSM_List
    fitsfile, prefix, nm, config = interpolate_locals("fitsfile prefix nm config")

    nm = int(nm)

    do_step = lambda step: step>=float(start) and step<=float(stop) 
    cellsize = None
    npix = None
    
    config = config or CFG
    # Use parameter file settings
    if config:
        with open(config) as params_std:
            params = json.load(params_std)
    
        # Remove empty strings and coments, and convert unicode characters to strings
        for key in params.keys():
            if params[key] == "" or key in ["#"]:
                del params[key]
            elif isinstance(params[key],unicode):
                params[key] = str(params[key])

        global SYNTHESIS, SCAN, DIRECTION, OBSERVATORY, ANTENNAS, DEFAULT_IMAGING_SETTINGS, SEFD

        prefix = params["sim_id"] or "hi-inator"
        
        SYNTHESIS = params["synthesis"]
        SCAN = params["scanlength"]
        DTIME = params["dtime"]
        DIRECTION = params["direction"]
        DEFAULT_IMAGING_SETTINGS = params["use_default_imaging_settings"]

        fitsfile = "%s/%s"%(INDIR, params["sky_model"])
        OBSERVATORY = params["observatory"].lower()

        if OBSERVATORY.startswith("vla") or OBSERVATORY.startswith("jvla"):
            OBSERVATORY = "vla"

        if OBSERVATORY.startswith("ska1mid") :
            OBSERVATORY = "meerkat"

        elif OBSERVATORY in ["kat7","kat-7","kat_7"]:
            OBSERVATORY = "kat-7"

        ANTENNAS = "./observatories/"+_ANTENNAS[params["observatory"].lower()]

        nm = params["split_cube"]
        addnoise = params["addnoise"]
        SEFD = params["sefd"]
        cellsize = params["cellsize"]
        npix = params["npix"]
        
        psf = params["psf"]
        dirty = params["keep_dirty_map"]
        clean = params["clean"]
       
        
    get_pw = lambda fitsname: abs(pyfits.open(fitsname)[0].header['CDELT1'])

    im.cellsize = "%farcsec"%(cellsize or get_pw(fitsfile)*3600.)
    im.npix = npix or 2048

    (ra, dec),(freq0, dfreq, nchan), naxis = fitsInfo(fitsfile)
    dfreq = abs(dfreq)
    chunks = nchan//nm
    
    mss = ['%s-%04d.MS'%(prefix,d) for d in range(nm)]

    if do_step(1):
        x.sh("fitstool.py --unstack=$prefix:freq:$chunks $fitsfile")
        lsms = [ prefix+"-%04d.fits"%d for d in range(nm) ]
        
        with open(LSMLIST,'w') as lsm_std:
            makedir(OUTDIR)
            lsm_std.write(",".join(lsms))
            
        MS_LSM_List = ["%s,%s"%(a, b) for a, b in zip(mss, lsms) ]

        x.sh('rm -f $CLEANS $DIRTYS $PSFS')

        scalenoise = 1
        if addnoise:
            scalenoise = math.sqrt(SCAN/float(SYNTHESIS)) if SCAN < SYNTHESIS else 1.0
        
        global RECENTRE
        RECENTRE = False

        def  make_and_sim(ms_lsm_comb="$MS_LSM", options={}, **kw):
            msname,lsmname = II(ms_lsm_comb).split(',')
            v.MS = msname
            v.LSM = lsmname
            adaptFITS(lsmname)
            _simms()
            
            simsky(addnoise=addnoise, scalenoise=scalenoise, **options)
            image(restore=clean, dirty=dirty, psf=psf)

            if clean:
                _add(im.RESTORED_IMAGE, CLEANS)
            if dirty:
                _add(im.DIRTY_IMAGE, DIRTYS)
            if psf:
                _add(im.PSF_IMAGE, PSFS)

        pper('MS_LSM', make_and_sim)
        
        if clean:
            images = get_list(CLEANS)
            restored_image = outname="%s/%s_restored.fits"%(OUTDIR, prefix)
            # combine images into a single image
            im.argo.combine_fits(images, outname=restored_image, ctype='FREQ', keep_old=False)
            info("Resulting restored image is at is: $restored_image")
            if RUN_SOURCE_FINDER:
                lsm.sofia_search(restored_image, threshold=4, do_reliability=True)

        if dirty:
            images = get_list(DIRTYS)
            dirty_image = outname="%s/%s_dirty.fits"%(OUTDIR, prefix)
            # combine images into a single image
            im.argo.combine_fits(images, outname=dirty_image, ctype='FREQ', keep_old=False)
            info("Resulting dirty image is at is: $dirty_image")

        if psf:
            images = get_list(PSFS)
            psf_image = outname="%s/%s_psf.fits"%(OUTDIR, prefix)
            # combine images into a single image
            im.argo.combine_fits(images, outname=psf_image, ctype='FREQ', keep_old=False)
            info("Resulting psf image is at is: $psf_image")

        info("Deleting temporary files")
        for item in images:
            x.sh("rm -f $item")

        x.sh('rm -f $CLEANS $DIRTYS $PSFS')


def image(msname='$MS', lsmname='$LSM', remove=None, **kw):
    """ imaging/cleaning function """
   
    global NSRC

    msname,lsmname = interpolate_locals('msname lsmname')
    v.MS = msname
    v.LSM = lsmname

    if DEFAULT_IMAGING_SETTINGS:
        # choose cellsize to be sixth of resolution
        im.npix = 2048
        im.weight = "briggs"
        im.robust = 2 # this equivalent to natural weighting

    ms.set_default_spectral_info()

    im.make_image(restored_image=RESTORED, 
            restore_lsm=False, mgain=0.85, fitbeam=True,**kw)

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


def extract_from_skymodel(lsmname='$LSM', addnoise=True, withnoise='$WITHNOISE', checkfirst=True,**kw):
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
                addnoise_fits(lsmname, withnoise,1e-9)

    lsm.sofia_search(withnoise, makeplot=True, **kw)

    return withnoise


def addnoise_fits(fitsname, outname, rms, correlated=False, clobber=True):
    """ add noise to fits file """

    im.argo.swap_stokes_freq(fitsname)
    hdu = pyfits.open(fitsname)
    hdr0 = hdu[0].header
    data = hdu[0].data
    noise =  np.random.normal(0, rms, data.shape) 

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

    
def _simms(msname='$MS', observatory='$OBSERVATORY', antennas='$ANTENNAS',
           synthesis=None, dtime=60, freq0='$FREQ0', dfreq='10MHz',
           nchan=2, direction='$DIRECTION',
           fromfits=True, fitsfile='$LSM',
           **kw):
    """ create simulated measurement set """

    msname, observatory, antennas, direction, freq0, fitsfile = \
        interpolate_locals('msname observatory antennas direction freq0 fitsfile')

    if MS_REDO is False and exists(msname):
        v.MS = msname
        return

    makedir(MSOUT)

    if fromfits:
        # I only need freq,dfreq,nchan from the bellow line
        (ra, dec), (freq0, dfreq, nchan), naxis = fitsInfo(fitsfile)
        freq0 += dfreq/2

        if SPWIDS > 1:
            bw = nchan*dfreq
            freq0 = [f0-dfreq/2 for f0 in np.arange(freq0,freq0+bw,bw/2)]
            _nchan = nchan/SPWIDS if nchan%SPWIDS==0 else (nchan +(SPWIDS-1)) /SPWIDS
            nchan = [_nchan]*SPWIDS
        dfreq = [dfreq]*SPWIDS

    synthesis = SCAN

    info('Making simulated MS: $msname. This may take a while')
    simms.create_empty_ms(msname, tel=observatory, direction=direction, 
                pos=antennas, nbands=SPWIDS, synthesis=synthesis, 
                dtime=dtime, freq0=freq0, dfreq=dfreq, nchan=nchan, **kw)

    v.MS = msname


def recentre_fits(fitsname, ra, dec):
    """ recentre FITS file """
    with pyfits.open(fitsname) as hdu:
        hdr = hdu[0].header
        hdr['CRVAL1'] = ra
        hdr['CRVAL2'] = dec
        hdu.writeto(fitsname, clobber=True)


def make_pure_lsm():
    makedir(DESTDIR)
    Model = galProps.Model(PURE_CAT, LSM)
    model = Model.load()
    model.writeto(PURE_CAT_LSM, overwrite=True)
    return model.nsrcs

def compute_vis_noise (sefd=None):
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
        if hdr['CTYPE3'].startswith('FREQ') is False: 
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

    for row0 in range(0, nrows, rowchunk):
        nr = min(rowchunk, nrows-row0)
        dshape[0] = nr
        data = noise*(numpy.random.randn(*dshape) + 1j*numpy.random.randn(*dshape)) * scale_noise

        if addToCol: 
            data+=colData[row0:(row0+nr)]
            info(" $addToCol + noise --> $column (rows $row0 to %d)"%(row0+nr-1))
        else : info("Adding noise to $column (rows $row0 to %d)"%(row0+nr-1))

        tab.putcol(column, data, row0, nr)
    tab.close() 


def adaptFITS(image):
    """ Try to re-structre FITS file so that it conforms to lwimager standard """
    
    hdr = pyfits.open(image)[0].header
    naxis = hdr["NAXIS"]
    
    # Figure if any axes have to be added be we proceed
    if naxis>=2 and naxis < 4:
        _freq = "--add-axis=freq:$FREQ0:DFREQ:Hz $image"
        _stokes = "--add-axis=stokes:1:1:1"
        if naxis == 2:
            info("FITS Image has 2 axes. Will add FREQ and STOKES axes. We need these predict visibilities")
            x.sh("fiitstool.py ${_stokes} ${_freq} $image")

        elif naxis==3:
            if hdr["CTYPE3"].lower().startswith("freq"):
                # Will also need to reorder if freq is 3rd axis
                x.sh("fitstool.py ${_stokes} --reorder=1,2,4,3 $image")

            elif hdr["CTYPE3"].lower().startswith("stokes"):
                x.sh("fitstool.py ${_freq} $image")

    info("You image is now LWIMAGER approved ;) ")
    


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
