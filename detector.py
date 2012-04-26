import sys
sys.path.append("utils")
import tools,pylab


class Detector:
    """ Area detector.
    """
    def __init__(self,parent=None):
        self._parent = parent
        # default values
        self.distance = 0.15
        self.pixelsize = 15E-06
        self.binning = 8
        self._Nx = 4096
        self._Ny = 4096
        # no saturation level defined
        self.saturationlevel = -1
        # init detector map
        self.detector_map = Detector_Map(x_gapsize_in_pixel=127,Nx=4096,Ny=4096)

    def get_effective_pixelsize(self):
        return self.pixelsize*self.binning         

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

class Detector_Map:

    def __init__(self,**kwargs):
        """
        Possible kwargs:
        
         - mask_array: array that defines the mask (0: masked out, 1: not masked out)

         - Nx: horizontal dimension of the mask (inside the mask, without counting any gaps)
         - Ny: vertical dimension of the mask (inside the mask, without counting any gaps)

         - x_gapwidth_in_pixel: horizontal gap is generated with given width (in unit unbinned pixel)
         - y_gapwidth_in_pixel: vertical gap is generated with given width (in unit unbinned pixel)

         - cx: x-coordinate of position (primary beam)
         - cy: y-coordinate of position (primary beam)

        """
        # init mask array
        if 'mask_array' in kwargs: 
            self.mask = kwargs['mask_array']
        elif 'Nx' in kwargs and 'Ny' in kwargs:
            Nx = kwargs['Nx']; Ny = kwargs['Ny']
            if 'x_gapwidth_in_pixel' in kwargs: Ny += kwargs['x_gapwidth_in_pixel']
            if 'y_gapwidth_in_pixel' in kwargs: Nx += kwargs['y_gapwidth_in_pixel']
            self.mask = pylab.ones(shape=(Ny,Nx))
            
        # set pixels in gap to zero
        if 'x_gapwidth_in_pixel' in kwargs: self.mask[pylab.ceil(self.cy)-kwargs['x_gapwidth_in_pixel']/2:pylab.ceil(cy)-kwargs['x_gapwidth_in_pixel']/2+kwargs['x_gapwidth_in_pixel'],:]
        if 'y_gapwidth_in_pixel' in kwargs: self.mask[:,pylab.ceil(self.cx)-kwargs['y_gapwidth_in_pixel']/2:pylab.ceil(cx)-kwargs['y_gapwidth_in_pixel']/2+kwargs['y_gapwidth_in_pixel']]

        # set center position
        if 'cx' in kwargs: self.cx = cx
        else: self.cx = (self.mask.shape[1]-1)/2.
        if 'cy' in kwargs: self.cy = cy
        else: self.cy = (self.mask.shape[0]-1)/2.
        
    def get_cx(binning=1): return (self.cx-(.5)/(1.*binning)

    def get_cy(binning=1): return self.cy/(1.*binning)
