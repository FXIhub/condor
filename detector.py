import sys
sys.path.append("utils")
import tools,pylab


class Detector:
    """ Area detector.
    """
    def __init__(self,parent=None):
        self.distance = 0.15
        self.pixelsize = 15E-06
        self.binning = 8
        self._Nx = 4096
        self._Ny = 4096
        self.saturationlevel = -1
        self._parent = parent
        gappixel = 127
        self.set_mask(gappixel*self.pixelsize,'x')

    def set_center(cx=None,cy=None,scaling='pixel'):
        if not cx: (self.mask.shape[1]-1)/2.0 
        else: self._cx = cx
        if not cy: (self.mask.shape[0]-1)/2.0 
        else: self._cy = cy

    def get_effective_pixelsize(self):
        return self.pixelsize*self.binning         

    def set_mask(self,gapsize,gaporientation,cx=None):
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
        if self.cx == None: cx = round(self.mask.shape[1]/2.0)
        else: cx = self.cx
        if self.cy == None: cy = round(self.mask.shape[0]/2.0)
        else: cy = self.cy
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

class Mask:

    def __init__(self,**kwargs):
        # generate mask array
        if 'mask' in kwargs: 
            self.mask = kwargs['mask']
        elif 'Nx' in kwargs and 'Ny' in kwargs:
            Nx = kwargs['Nx']; Ny = kwargs['Ny']
            if 'gapwidth_in_pixel' in kwargs:
                if kwargs['gap_orientation'] == 'x': Ny += kwargs['gapwidth_in_pixel']
                elif kwargs['gap_orientation'] == 'y': Nx += kwargs['gapwidth_in_pixel']
                else: print "Warning: No valid gap orientation given while initializing mask."
            self.mask = pylab.ones(shape=(Ny,Nx))
            
        # set center position
        if 'cx' in kwargs: self.cx = cx
        else: self.cx = (Nx-1)/2.
        if 'cy' in kwargs: self.cy = cy
        else: self.cy = (Ny-1)/2.

        if 'gapwidth_in_pixel' in kwargs:
            # set pixels in gap to zero
            if kwargs['gap_orientation'] == 'x':
                self.mask[pylab.ceil(self.cy)-kwargs['gapwidth_in_pixel']/2:pylab.ceil(cy)-kwargs['gapwidth_in_pixel']/2+kwargs['gapwidth_in_pixel'],:]
                self.mask[:,pylab.ceil(self.cx)-kwargs['gapwidth_in_pixel']/2:pylab.ceil(cx)-kwargs['gapwidth_in_pixel']/2+kwargs['gapwidth_in_pixel']]
        
