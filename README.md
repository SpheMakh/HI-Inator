# HI-Inator
Radio interefometry simulator/imager tailored for HI sky models.

# Download
```
git clone https://github.com/SpheMakh/HI-Inator
```

# Available Arrays
Request other Arrays via the [issues](https://github.com/SpheMakh/HI-Inator/issues) service.   
Specify the array you want to use by its key in config file (*/input/parameters.json*).  

| Key | Array |    
| ------|-----|  
|meerkat|MeerKAT|  
|kat-7|KAT-7|  
|jvlaa|JVLA A Config|  
|jvlab|JVLA B Config|  
|jvlac|JVLA C Config|    
|jvlad|JVLA D Config|  
|wsrt|WSRT|  
|ska1mid254|SKA1MID 254 dishes|  
|ska1mid197|SKA1MID 197 dishes|  


# Running the simulator
1. Place your input skymodel in HI-Inator/input  
2. Modify the *input/parameters.json* config file to suite your needs

## Without Docker (Not recommended)
To run HI-Inator without docker you will need to install MeqTrees and all its related software. This is much easier these days (Thanks to Gijs Molenaar), see [radio-astro ppa](https://launchpad.net/~radio-astro/+archive/ubuntu/main). These are the packages that will need (I may niss some of them):

1. Meqtrees Timba
2. Meqtrees Cattery, Owlcat, Purr  
3. Pyxis
4. LWIMAGER
5. WSClean (Optional but very usefull)
6. Casacore, casacore data
7. CASAPY

Once you have all these installed: you can simply run the pipeline as follows. 
```
cd HI-Inator/src
pyxis CFG=../input/parameters.json OUTDIR=../output azishe
```

Then you simulation reseults--images and visibilities (Measurement set)-- should be in HI-Inator/output

## With Docker (Recommended)
First make sure you have the latest docker (>= 1.3) installed (not the default Ubuntu docker.io package).

https://docs.docker.com/installation/ubuntulinux/#ubuntu-trusty-1404-lts-64-bit

Once you have docker setup, run: 

1. Download casapy `$ make download`  
2. Build container `$ make build`  
3. Run simulation `$ make run`  
