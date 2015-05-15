#!/usr/bin/env python
## Reads and manipulates Ed's galaxy property catalogs
## Requires pyfits, astLib
## sphe spehmakh@gmail.com

# import python essentials
import os
import sys
import time
import numpy as np
import math
import tempfile
import subprocess

# non-standard packages
import pyfits
from astLib.astWCS import WCS
import Tigger

# some useful constants
PI = math.pi
FWHM = math.sqrt(math.log(256))


MAP = {'Galaxy':'ID',
'ALOG10(MHI)':'mhi',
'lwidth':'lwidth',
'incl':'incl',
'z_obs':'z_obs',
'dist [Mpc]':'dist',
'Sint [Jy km/s]':'int_flux',
'RA-DEC pixel position in master cube':'radec_pix',
'Scaled cubelet dimensions':'scd',
'dV [km/s]':'dV',
'chan1, chan2':'chans',
}

def mysplit(string,delimiter=None):
    if delimiter:
        return [item.strip() for item in string.split(delimiter)]
    else:
        return [item.strip() for item in string.split()]

# Communication functions
def info(string):
    t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    print "%s ##INFO: %s"%(t,string)
def warn(string):
    t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    print "%s ##WARNING: %s"%(t,string)
def abort(string):
    t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    raise SystemExit("%s ##ABORTING: %s"%(t,string))


# run things on the command line
def _run(command,options):
    cmd = " ".join([command]+options)
    info('running: %s'%cmd)
    process = subprocess.Popen(cmd,
                  stderr=subprocess.PIPE if not isinstance(sys.stderr,file) else sys.stderr,
                  stdout=subprocess.PIPE if not isinstance(sys.stdout,file) else sys.stdout,
                  shell=True)
    if process.stdout or process.stderr:
        out,err = process.comunicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        out = None
    else:
        process.wait()
    if process.returncode:
            abort('%s: returns errr code %d'%(command,process.returncode))


class Model(object):
    """
    Class that reads and manipulates Ed's galaxy property files
    """
    def __init__(self,textname,fitsname,sources=None,ra0=0,dec0=0):
        """
        fitsname : Name of FITS cube where galaxies are simulated
        """
        if sources is None:
            sources = []
        self.textname = textname
        self.fitsname = fitsname
        self.ra0 = ra0
        self.dec0 = dec0
        self.sources = sources
        self.nsrcs = len(sources)

    def load(self,textname=None,fitsname=None,append=False):
        """ load galaxy properties from text file """

        if textname is None:
            textname = self.textname
        if fitsname is None: 
            fitsname = self.fitsname
        if append:
            if isinstance(self.textname,list):
                self.textname.append(textname)
            else:
                self.textname = [self.textname,textname]
            sources = self.sources
        else:
            sources = []
        self.fitsname = fitsname

        nsrc = 0
        info('loading galaxy properties from %s'%textname)
        info('Corresponding FITS cube is %s'%fitsname)

        std = open(textname)
        wcs = WCS(pyfits.open(fitsname)[0].header,mode='pyfits')    
        nx,ny = pyfits.open(fitsname)[0].data.shape[-2:]
        self.ra0,self.dec0 = wcs.pix2wcs(nx/2.,ny/2.)

        for line in std.readlines():
            if line[0] not in ['\n','#']:
                key,val = mysplit(line,':')
                if key == 'Galaxy':
                    nsrc +=1
                    src = Source()
                    sources.append(src)
                if key.startswith('Scaled'):
                    val = mysplit(val)
                elif key.startswith('chan'):
                    chan1,chan2 = mysplit(val)
                    src.addAttribute('chan1',chan1)
                    src.addAttribute('chan2',chan2)
                elif key.startswith('RA-DEC'):
                    ra_pix,dec_pix = map(float,mysplit(val))
                    src.addAttribute('ra_pix',ra_pix)
                    src.addAttribute('dec_pix',dec_pix)
                    ra_deg,dec_deg = wcs.pix2wcs(ra_pix,dec_pix)
                    src.addAttribute('ra_deg',ra_deg)
                    src.addAttribute('dec_deg',dec_deg)
                else:
                    src.addAttribute(MAP[key],val)

        self.sources = sources
        self.nsrcs = nsrc
        info('Loaded %d sources'%nsrc)
        return self

    def writeto(self,filename,overwrite=False):
        
        if os.path.exists(filename) and overwrite is False:
            abort('%s already exists. Set overwrite to True to overwrite')
        tf = tempfile.NamedTemporaryFile(suffix='.txt')
        tf.write('#format:name ra_d dec_d emaj_s emin_s pa_d i\n')
        for i,src in enumerate(self.sources):
            tf.write('%d %.8g %.8g 0 0 0 %.4g\n'%(src.ID,src.ra_deg,src.dec_deg,src.int_flux/src.dV))
        tf.flush()

        options = [tf.name,filename]
        if overwrite:
            options += ['-f']
        _run('tigger-convert',options=options)
        tf.close()

        # add non-standard attributes as tags
        model = Tigger.load(filename)
        for src0,src1 in zip(sorted(self.sources,key=lambda src:src.ID),
                             sorted(model.sources,key=lambda src: int(src.name))):
            for attribute in 'mhi dV dist incl z_obs scd chan1 chan2'.split():
                src1.setAttribute(attribute,getattr(src0,attribute))
        model.save(filename)
        # Rename sources according to the COPART
        # _run('tigger-convert',[filename,'--rename','-f','--min-extent','0','--cluster-dist','10'])
            

class Source(object):

    def __init__(self):
        """
        ID : Galaxy ID
        int_flux : Integrated flux (stokes I)
        mhi: HI mass
        lwidth : Line width
        incl : Inclination angle
        obs_z : Observed redshift
        dist : Distance in Mpc
        ra_pix : RA in pixels
        dec_pix : DEC in pixels
        ra_deg : RA in degrees
        dec_deg : DEC in degrees
        scd : Scaled cubelet dimension
        dV : Velocity
        chan1 : Initial channel
        chan2 : Fianl channel
        """

        self.ID = None
        self.int_flux = None
        self.mhi = None
        self.lwidth = None
        self.incl = None
        self.obs_z = None
        self.dist = None
        self.ra_pix = None
        self.dec_pix = None
        self.ra_deg = None
        self.dec_deg = None
        self.scd = None
        self.dV = None
        self.chan1 = None
        self.chan2 = None

    def addAttribute(self,key,val):
        try:
            val = float(val)
            if val%1 == 0:
                val = int(val)
        except ValueError:
            val = val
        setattr(self,key,val)
