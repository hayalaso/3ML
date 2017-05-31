#Reading a Galprop Map
from threeML.models.spatialmodel import SpatialModel
from threeML.models.Parameter import Parameter,SpatialParameter
import pyfits as pf
import numpy as np
from astropy import coordinates as coord
from astropy import units as u

def return_lonlat(RA,Dec,sky="CEL"):
    #sky tells in what coordinate system is the position
    #Map is in Galactic coords, need to convert from Celestial
    if sky=="CEL":
        c = coord.SkyCoord(ra=RA, dec=Dec, unit=(u.degree, u.degree),
                          frame='icrs')
        lon = c.galactic.l.degree
        lat = c.galactic.b.degree
    else:
        lon, lat = RA, Dec

    return lon, lat

def coordToPixel(l,b,filename):
    m=pf.open(filename)
    stl = m[0].header['CRVAL1']
    dl = m[0].header['CDELT1']
    if 'CRPIX1' in m[0].header:
        pix0l = m[0].header['CRPIX1']
    else:
        pix0l = 0

    stb = m[0].header['CRVAL2']
    db = m[0].header['CDELT2']
    if 'CRPIX2' in m[0].header:
        pix0b = m[0].header['CRPIX2']
    else:
        pix0b = 0

    il = int((l-stl)/dl + pix0l)
    ib = int((b-stb)/db + pix0b)

    return np.array([ib, il])

def fp_linear_interp(pos,dat,filename):
        #4-point linear interpolation to determine the values at the requested coords
        px_0=pos[-1].astype('int')
        py_0=pos[-2].astype('int')
        px_arr=np.array([px_0,px_0,px_0+1,px_0+1])
        py_arr=np.array([py_0,py_0+1,py_0,py_0+1])
        b=pos[-1]-px_0
        a=1.-b
        d=pos[-2]-py_0
        c=1.-d
        fp_pos=[py_arr,px_arr]
        try:
            vals = dat[pos[-3]][fp_pos]
        except ValueError:
            print "The Galprop map in file {} is not defined at the sky coordinates requested by the user".format(filename)


        vals=vals*np.array([a*c,a*d,b*c,b*d]) #total area is 1, so factors are already normalized
        vals=np.sum(vals,axis=0)

        return vals

class GalPropMap(SpatialModel):
    """Read a Galprop map and get expectations for diffuse emission """
    def setup(self,file):
        #assumes that file is a fits file with data in hdu 0
        #flux map is in units of (MeV s sr cm^2)^-1
        self.filename=file
        self.w=pf.open(self.filename)
        self.data=np.nan_to_num(pf.getdata(file,0))
        self.cube=self.data[0] #Galprop map usually has dimension [1,35,360,720]
        #energies are in log(E) E is in MeV
        self.energies=np.array(self.w[0].header['CRVAL3']+i*self.w[0].header['CDELT3'] for i in range(self.w[0].header['NAXIS3']))
        self.functionName        = "GalPropSpatialMap"
        self.parameters          = collections.OrderedDict()
        self.parameters['Norm']  = Parameter('Norm',1.,0.,1.e5,0.1,fixed=True,nuisance=False,dataset=None)
        self.ncalls              = 0

    def __call__(self,RA,Dec,energy):
        self.ncalls += 1
        Norm = self.parameters['Norm'].value
        #Map is in Galactic coords need to convert from Celestial
        lon, lat = return_lonlat(RA,Dec)
        #Get the pixel position
        px, py = coordToPixel(lon,lat,self.filename)
        #determine values at requested energie by PL interp on the nearest available in the model
        E_0=np.max(np.where(energies<=np.log10(energy)),axis=1)[0]
        pos_0=[E_0,py,px]
        vals_0=fp_linear_interp(pos_0,self.cube,self.filename)
        pos_1=[E_0+1,py,px]
        vals_1=fp_linear_interp(pos_1,self.cube,self.filename)
        gamma=-(np.log10(vals_1)-np.log10(vals_0))/self.w[0].header['CDELT3']
        vals=vals_0-gamma*(np.log(energy)-logE0)
        vals=np.pow(10,vals)
        return Norm*vals
