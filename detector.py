import constants as phy
import pylab
import tools

class Detector:
    """ Area detector.
    """
    def __init__(self,parent=None):
        self.distance = 0.15
        self.pixelsize = 15E-06
        self.binning = 8
        self.Nx = 4096
        self.Ny = 4096
        self.saturationlevel = -1
        self._parent = parent
        gappixel = 127
        self.set_mask(gappixel*self.pixelsize,'x')

    def get_effective_pixelsize(self):
        return self.pixelsize*self.binning         

    def set_mask(self,gapsize,gaporientation):
        """ Masks out regions within the gap, creates mask of the size of the later generated image. "0" denotes masked out pixels, "1" active pixels.
        """
        Nx_eff = round(self.Nx/self.binning)
        Ny_eff = round(self.Ny/self.binning)
        gappixel = round(gapsize/self.get_effective_pixelsize())
        if gaporientation == "x":
            Ny_eff += gappixel
            self.mask = pylab.ones(shape=(Ny_eff,Nx_eff))
            self.mask[round(Ny_eff/2.0-0.5-gappixel/2.0):round(Ny_eff/2.0-0.5-gappixel/2.0)+gappixel,:] = 0
        elif gaporientation == "y":
            Nx_eff += gappixel
            self.mask = pylab.ones(shape=(Ny_eff,Nx_eff))
            self.mask[:,round(Nx_eff/2.0-0.5-gappixel/2.0):round(Nx_eff/2.0-0.5-gappixel/2.0)+gappixel] = 0


    # Functions that convert detector-coordinates:
    def _get_q_from_r(self,r):
        return 4*pylab.pi/self._parent.source.photon.get_wavelength()*pylab.sin(pylab.arctan(r/self.distance)/2.0)
        #return 2.0*pylab.pi/self._parent.source.photon.get_wavelength()*r/self.distance
    
    def _get_r_from_q(self,q):
        return self.distance*pylab.tan(pylab.arcsin(q*self._parent.source.photon.get_wavelength()/4.0/pylab.pi)*2.0)
        #return self.distance*q*self._parent.source.photon.get_wavelength()/2.0/pylab.pi
    
    def _get_i_from_q(self,q):
        return int(round(self._get_r_from_q(self,q)/self.get_effective_pixelsize()))
   
    def _get_q_from_i(self,i):
        return self._get_q_from_r(self,i*self.get_effective_pixelsize())

    def write_qmap(self,wavelength=None):
        if not wavelength:
            wavelength = self._parent.source.photon.get_wavelength()
        qmap = pylab.ones_like(self.mask)
        X,Y = pylab.meshgrid(pylab.arange(0,self.mask.shape[1]),pylab.arange(0,self.mask.shape[0]))
        cx = round(self.mask.shape[1]/2.0)
        cy = round(self.mask.shape[0]/2.0)
        X -= cx
        Y -= cy
        R = pylab.sqrt(X**2+Y**2)
        qmap = self._get_q_from_r(R*self.get_effective_pixelsize())
        return qmap

    # Crystallographic resolution (full-period resolution at the closest edge)
    def get_max_crystallographic_resolution(self):
        min_detector_center_edge_distance = min([self.mask.shape[1]-self.cx,self.cx,self.mask.shape[0]-self.cy,self.cy])*self.pixelsize
        wavelength = self._parent.source.photon.get_wavelength()
        detector_distance = self.distance
        return tools.get_max_crystallographic_resolution(wavelength,min_detector_center_edge_distance,detector_distance)
