## HI Simulations pipeline
## github.com:SpheMakh/Ceiling-KAT_Apps [branch hikat]
## sphe sphemakh@gmail.com

# import some python essentials
import os
import sys
import math
import numpy as np
import tempfile
#import mqt
import json
import scipy.ndimage as nd
import subprocess

# some useful constants
PI = math.pi
FWHM = math.sqrt(math.log(256))

# import pyxis essentials
import Pyxis
from Pyxis.ModSupport import *
import im
import ms
import lsm


# import non standard package
import pyfits
import Tigger
import galProps
import pyrap.measures 
dm = pyrap.measures.measures()
import pyrap.quanta as dq


import im.lwimager 


def simsky(gain_err=None, pbeam=False,
           pointing_err=None, addnoise=True, sefd=551,
           noise=0, scalenoise=1, column='$COLUMN',
           freq_chunks='$FREQ_CHUNKS',
           **kw):
    """ simulate sky into MS """

    column, freq_chunks = interpolate_locals('column freq_chunks')

    if RECENTRE:
        ra, dec = map(np.rad2deg, [ dm.direction(*DIRECTION.split(","))[m]["value"] for m in "m0","m1"])
        recentre_fits(LSM, ra=ra, dec=dec)

    ms.set_default_spectral_info()
    im.IMAGE_CHANNELIZE = 1

    im.lwimager.predict_vis(image=LSM, column="MODEL_DATA", padding=1.5, wprojplanes=128)
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


def try_again(cmd, exception, tries=3):
    for i in range(tries):
        try:
            cmd()
            return
        except exception:
            warn(" Attempt $i of running command $cmd failed. %d attemps left"%(tries-i-1) )
    raise exception


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
    keep_ms = False
    

    # Initialise CASAPY if running on container
    proc = subprocess.Popen(["sh","inDocker.sh"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.stderr.read():
        info("Cannot tell whether I'm running on a container or not. Will initialise CASA just as precaution")
        x.sh("casapy --help --log2term --nologger --nogui --help -c quit")

    elif proc.stdout.read().strip() == "True":
        info("Initialising CASA on the container")
        x.sh("casapy --help --log2term --nologger --nogui --help -c quit")
    
    config = config or CFG
    source_finder = {}
    component_model = False
    # Use parameter file settings
    if config:
        with open(config) as params_std:
            jparams = json.load(params_std)
    
        # create a new dict from the json file
        params = {}
        for key, val in jparams.iteritems():
            # Make sure all keys are strings
            _key = str(key)
            
            # ignore empty strings and comments
            if val=="" or _key=="#":
                pass
            #convert unicode values to strings
            elif isinstance(val,unicode):
                params[_key] = str(val)
            else:
                params[_key] = val


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
        ncpu = params.get("ncpu", nm)
        addnoise = params["addnoise"]
        SEFD = params["sefd"]
        cellsize = params["cellsize"]
        npix = params["npix"]
        
        psf = params["psf"]
        dirty = params["keep_dirty_map"]
        clean = params["clean"]
        keep_ms = params["keep_ms"]
        component_model = params.get("component_model", False)
    
        # Source finder business
        global RUN_SOURCE_FINDER
        RUN_SOURCE_FINDER = params["run_source_finder"]
        source_finder["threshold"] = params["sf_thresh"]

        if params["sf_conf"]:
            sf_conf = params["sf_conf"]
            look_here = [sf_conf] + map(II, "input/${input:FILE} ${input:FILE}".split() )

            found_path = False
            for _path in look_here:
                if exists(_path):
                    source_finder["sofia_conf"] = _path
                    break
            if not found_path:
                warn("Cannot find your custom Sofia conf. Will use system default")
            
        ncores(ncpu or nm)

    get_pw = lambda fitsname: abs(pyfits.open(fitsname)[0].header['CDELT1'])

    im.cellsize = "%farcsec"%(cellsize or get_pw(fitsfile)*3600.)
    im.npix = npix or 2048

    (ra, dec),(freq0, dfreq, nchan), naxis = fitsInfo(fitsfile)
    dfreq = abs(dfreq)
    chunks = nchan//nm
        
    v.MS_List = mss = ['%s/%s-%04d.MS'%(MSOUT,prefix,d) for d in range(nm)]

    if do_step(1):
        x.sh("fitstool.py --unstack=$prefix:freq:$chunks $fitsfile")
        lsms = [ prefix+"-%04d.fits"%d for d in range(nm) ]
        
        with open(LSMLIST,'w') as lsm_std:
            makedir(OUTDIR)
            lsm_std.write(",".join(lsms))
            
        MS_LSM_List = ["%s,%s"%(a, b) for a, b in zip(mss, lsms) ]

        x.sh('rm -f $CLEANS $DIRTYS $PSFS $MODELS $RESIDUALS')

        scalenoise = 1
        if addnoise:
            scalenoise = math.sqrt(SCAN/float(SYNTHESIS)) if SCAN < SYNTHESIS else 1.0
        
        global RECENTRE
        RECENTRE = True
    
        def  make_and_sim(ms_lsm_comb="$MS_LSM", options={}, **kw):
            msname,lsmname = II(ms_lsm_comb).split(',')
            v.MS = msname
            v.LSM = lsmname
            adaptFITS(lsmname)
            _simms()
            
            simsky(addnoise=addnoise, scalenoise=scalenoise, **options)

#            if component_model:
#                mqt_opts = {"sim_mode":"add to MS", 
#                             "tiggerlsm.filename": component_model,
#                             "ms_sel.output_column": "CORRECTED_DATA"}
#
#                mqt.msrun(TURBO_SIM, job='_tdl_job_1_simulate_MS', section="sim", options=mqt_opts)

            image(restore=clean, dirty=dirty, psf=psf)

            if clean:
                _add(im.MODEL_IMAGE, MODELS)
                _add(im.RESIDUAL_IMAGE, RESIDUALS)
                _add(im.RESTORED_IMAGE, CLEANS)
            if dirty:
                _add(im.DIRTY_IMAGE, DIRTYS)
            if psf:
                _add(im.PSF_IMAGE, PSFS)
        
        pper('MS_LSM', make_and_sim)
        v.MS = "%s/%s.MS"%(OUTDIR, prefix)
        im.argo.icasa("concat", vis=mss, concatvis=MS, timesort=True)
        # clean up

        for msname in mss:
            x.sh("rm -fr $msname")

        #ms.virtconcat(MS, thorough=True)
        # compress and/or delete MS/s
        if keep_ms and config:
            x.sh("tar -czf ${MS}.tar.gz $MS && rm -fr $MS")
        else:
            x.sh("rm -fr $MS")

        # combine images into a single image
        def _combine(imlist, label):
            images = get_list(imlist)
            image_name = outname="%s/%s_%s.fits"%(OUTDIR, prefix, label)
            im.argo.combine_fits(images, outname=image_name, ctype='FREQ', keep_old=False)
            return image_name

        
        if clean:
            names = []
            for imlist, label in zip([MODELS, RESIDUALS, CLEANS],["model", "residual", "restored"]):
               names.append( _combine(imlist, label) )
            info("Resulting clean,model and residual images are : $restored_image, %s %s %s"%tuple(names))

            if RUN_SOURCE_FINDER:
                lsm.sofia_search(fitsname=names[-1], do_reliability=False, **source_finder)

        if dirty:
            info("Resulting dirty image is at: %s"%_combine(DIRTYS, "dirty"))

        if psf:
            info("Resulting psf image is at: %s"%_combine(PSFS, "psf"))

        info("Deleting temporary files")
        x.sh("rm -f ${OUTDIR>/}*wsclean*-first-residual.fits") # Not sure why wsclean produces these
        x.sh("rm -f ${OUTDIR>/}*wsclean*-MFS*.fits") # Remove MFS images
        x.sh('rm -f $MODELS $RESIDUALS $CLEANS $DIRTYS $PSFS')
        info("DONE!")


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

    im.make_image(restore_lsm=False, mgain=0.85, fitbeam=True,**kw)

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
    im.argo.make_image(psf=True,dirty=False,psf_image=tfname)
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
    makedir(MSOUT)

    msname, observatory, antennas, direction, freq0, fitsfile = \
        interpolate_locals('msname observatory antennas direction freq0 fitsfile')
    
    if MS_REDO is False and exists(msname):
        v.MS = msname
        return


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
    ms.create_empty_ms(msname, tel=observatory, direction=direction, 
                pos=antennas, nbands=SPWIDS, synthesis=synthesis, 
                dtime=dtime, freq0=freq0, dfreq=dfreq, nchan=nchan, **kw)

    v.MS = msname


def recentre_fits(fitsname, ra, dec):
    """ recentre FITS file """
    with pyfits.open(fitsname) as hdu:
        hdu[0].header['CRVAL1'] = ra
        hdu[0].header['CRVAL2'] = dec
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
    noise = sefd/math.sqrt(abs(2*bw*dt))
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
                x.sh("fitstool.py ${_stokes} $image && fitstool.py --reorder=1,2,4,3 $image")

            elif hdr["CTYPE3"].lower().startswith("stokes"):
                x.sh("fitstool.py ${_freq} $image")


    with pyfits.open(image) as hdu:
        hdr = hdu[0].header
        freq_ind = filter(lambda ind: hdr["CTYPE%d"%ind].startswith("FREQ"), range(1,5) )[0]

        if hdr["CDELT%d"%freq_ind] <0:
            dfreq = abs(hdr["CDELT%d"%freq_ind])
            nchan = hdr["NAXIS%d"%freq_ind]
            freq0 = hdr["CRVAL%d"%freq_ind]
            hdu[0].header["CDELT%d"%freq_ind] = dfreq
            hdu[0].header["CRVAL%d"%freq_ind] = freq0 - nchan*dfreq
            hdu.writeto(image, clobber=True)
        
            
    info("You image is now OMS approved ;) ")

    

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
