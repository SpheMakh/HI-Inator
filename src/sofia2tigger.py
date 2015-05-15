#!/usr/bin/env python
## converts a sofia catalog to a tigger lsm. Requires: Tigger, tigger-convert
## sphe sphemakh@gmail.com

import Tigger
import os
import sys
import subprocess
import tempfile
from argparse import ArgumentParser
import numpy as np

mapping = {'ID':'name',
'ID_old':'ID_old',
'Xg':'Xg',
'Yg':'Yg',
'Zg':'Zg',
'Xm':'Xm',
'Ym':'Ym',
'Zm':'Zm',
'Xmin':'Xmin',
'Xmax':'Xmax',
'Ymin':'Ymin',
'Ymax':'Ymax',
'Zmin':'Zmin',
'Zmax':'Zmax',
'NRvox':'NRvox',
'SNRmin':'SNRmin',
'SNRmax':'SNRmax',
'SNRsum':'SNRsum',
'NRpos':'NRpos',
'NRneg':'NRneg',
'Rel':'Rel',
'ELL_MAJ':'emaj_s',
'ELL_MIN':'emin_s',
'ELL_PA':'pa_d',
'F_PEAK':'F_PEAK',
'F_TOT':'i',
'F_Wm50':'F_Wm50',
'RMS_CUBE':'rms_cube',
'W20':'W20',
'W50':'W50',
'Wm50':'Wm50',
'RAm':'ra_d',
'DECm':'dec_d',
'FREQm':'freq0',
}

ADD_ATTRIBUTES = 'Xg,Yg,Zg,Rel,Xg,Yg,Zg,Xm,Ym,Zm,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,F_PEAK,SNRmin,SNRmax,SNRsum'

def _run(command,options):
    cmd = " ".join([command]+options)
    print 'running: %s'%cmd
    process = subprocess.Popen(cmd,
                  stderr=subprocess.PIPE if not isinstance(sys.stderr,file) else sys.stderr,
                  stdout=subprocess.PIPE if not isinstance(sys.stdout,file) else sys.stdout,
                  shell=True)
    if process.stdout or process.stderr:
        out,err = process.comunicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        out = None;
    else:
        process.wait()
    if process.returncode:
            print '%s: returns errr code %d'%(command,process.returncode)

if __name__=='__main__':

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
    parser = ArgumentParser(description='Convert SoFiA catalog to Tigger LSM (.lsm.html)')
    add = parser.add_argument
    add('-in','--input',
            help='input: Sofia Catalog')
    add('-out','--output',
            help='Tigger LSM output name (including extenstion)')
    add('-pw','--pixel-width',dest='pixel_width',type=float,
            help='Pixel width [in degrees]. Need this to '
                 'convert source dimensions to arcseconds')
    add('-a','--add',action='store_true',
            help='Add non-standard attributes (from sifia catalog) to the output catalog')
    add('-at','--attributes',default=ADD_ATTRIBUTES,
            help='Comma seperated attributes' 
                 ' to add to the output catalog. Enable --add. default: %s'%ADD_ATTRIBUTES)
    add('-me','--min-extent',dest='min_extent',default=0,type=float,
			help='see tigger-convert')
    add('-cd','--cluster-dist',dest='cluster_dist',default=0,type=float,
            help='see tigger-convert')
    args = parser.parse_args()

    filename = args.input
    output = args.output
    pixwidth = args.pixel_width

    std = open(filename,'rw')
    std.readline() # ignore first line
    names = std.readline().split()[1:]
    N = len(names)
    dtype = (' f64'*N).split()

    # First, lets convert source dimenstions from pixels to arcsecs
    sofia = np.genfromtxt(std,names=names,dtype=dtype)
    sofia['ELL_MAJ'] *= pixwidth*3600
    sofia['ELL_MIN'] *= pixwidth*3600
    # remove sources with zero flux [not sure why they are in the catalog]

    #TODO(sphe): Might be better to just create a fresh Tigger Model instead of using tigger-convert
    tf = tempfile.NamedTemporaryFile(suffix='.txt')
    tf.flush()
    tfname = tf.name
    np.savetxt(tfname,sofia)
    
    tigger_format = " ".join(["%s"%mapping[name] for name in names])
    _run('tigger-convert',['-t ASCII','--format=\"%s\"'%tigger_format,'-f','--min-extent',str(args.min_extent),'--cluster-dist',str(args.cluster_dist),tfname,output])
    tf.close()
    # add attributes
    if args.add:
        attributes = args.attributes.split(',')
        model = Tigger.load(output)
        for i,src in enumerate(sorted(model.sources,key=lambda src:src.name)):
            for attribute in attributes:
                src.setAttribute(attribute,sofia[attribute][i])
            if src.flux.I==0:
                model.sources.remove(src)
        model.save(output)
        # Rename sources according to COPART
        _run('tigger-convert',['--rename','-f',output,output])
