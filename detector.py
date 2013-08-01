import sys
sys.path.append("utils")
import pylab


class Detector:
    """
    A subclass of the input object.
    Defines area detector.
    """
    def __init__(self,**kwargs):
        """
        Function initializes Detector object.
        =====================================
        Arguments:
        Keyword arguments (if not given variable is set to default value):

        - distance: Distance between detector and interaction point. [0.741]
        - pixel_size: Edge length of square pixel (unbinned). [75E-06]
        - cx: Horizontal beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - cy: Vertical beam position in pixel. If argument is \'None\' or not given center is set to the middle. [None]
        - binning: Number of binned pixels (binning x binning). [1]
        - saturation_level: Detector saturation level in photons. \'saturation_level=None\' means that no saturation level is defined. [None]
        - parent: Input object that includes detector object. This variable is optional. [None]

        EITHER (default):
        - Nx: Number of pixels in horizontal direction not including a potential gap. [1024]
        - Ny: Number of pixels in vertical direction not including a potential gap. [1024]
        - x_gap_size_in_pixel: Width of central horizontal gap in pixel. [0]
        - y_gap_size_in_pixel: Width of central vertical gap in pixel. [0]
        
        OR:
        - mask: 2d mask array. In the array 1 stands for valid pixel and 0 for unvalid pixel (masked out).

        """

        self._parent = kwargs.get('parent',None)
        self.distance = kwargs.get('distance',0.15)
        self.pixel_size = kwargs.get('pixel_size',75.E-6)
        gx = kwargs.get('x_gap_size_in_pixel',0)
        gy = kwargs.get('y_gap_size_in_pixel',0)
        hd = kwargs.get('hole_diameter_in_pixel',0)
        self.binning = kwargs.get('binning',1)
        if 'mask' in kwargs:
            self.init_mask(mask=kwargs['mask'])
            self.cx = kwargs.get('cx',(self.mask.shape[1])/(2.*self.binning))
            self.cy = kwargs.get('cy',(self.mask.shape[0])/(2.*self.binning))   
        else:
            self.Nx = kwargs.get('Nx',1024)
            self.Ny = kwargs.get('Ny',1024)
            self.cx = kwargs.get('cx',(self.Nx+gy-1)/2.)
            self.cy = kwargs.get('cy',(self.Ny+gx-1)/2.)   
            self.init_mask(Nx=self.Nx,Ny=self.Ny,x_gap_size_in_pixel=gx,y_gap_size_in_pixel=gy,hole_diameter_in_pixel=hd,binning=self.binning)
        self.saturation_level = kwargs.get('saturation_level',None)

    def init_mask(self,**kwargs):
        """        
        Function initializes the detector mask.
        =======================================

        Arguments:
        
        Keyword arguments (if not given variable is set to default value):
        
        EITHER
        - mask: array that defines the mask (0: masked out, 1: not masked out)
        
        OR:
        - Nx: horizontal dimension of the mask (inside the mask, without counting any gaps)
        - Ny: vertical dimension of the mask (inside the mask, without counting any gaps)
        - x_gap_size_in_pixel: horizontal gap is generated with given width (in unit unbinned pixel)
        - y_gap_size_in_pixel: vertical gap is generated with given width (in unit unbinned pixel)
        - hole_diameter_in_pixel: holeis generated with given diameter (in unit unbinned pixel)

        """

        # init mask array
        if 'binning' in kwargs:
            self.binning = kwargs['binning']
        else:
            self.binning = 1
        if 'mask' in kwargs: 
            self.mask = kwargs['mask']
        elif 'Nx' in kwargs and 'Ny' in kwargs:
            Nx = kwargs['Nx']/self.binning; Ny = kwargs['Ny']/self.binning
            if 'x_gap_size_in_pixel' in kwargs:
                x_gap = kwargs['x_gap_size_in_pixel']/self.binning
                Ny += x_gap
            if 'y_gap_size_in_pixel' in kwargs: 
                y_gap = kwargs['y_gap_size_in_pixel']/self.binning
                Nx += y_gap
            self.mask = pylab.ones(shape=(Ny,Nx))
        else:
            print "ERROR: Either \'mask_array\' or \'Nx\' and \'Ny\' have to be specified."
            return 
            
        # set pixels in gap to zero
        if 'x_gap_size_in_pixel' in kwargs:
            if kwargs['x_gap_size_in_pixel'] != 0:
                cy = pylab.ceil((self.mask.shape[0]-1)/2.)
                self.mask[cy-kwargs['x_gap_size_in_pixel']/2/self.binning:cy-kwargs['x_gap_size_in_pixel']/2/self.binning+kwargs['x_gap_size_in_pixel']/self.binning,:] = 0
        if 'y_gap_size_in_pixel' in kwargs:
            if kwargs['y_gap_size_in_pixel'] != 0:
                cx = pylab.ceil((self.mask.shape[1]-1)/2.)
                self.mask[:,cx-kwargs['y_gap_size_in_pixel']/2/self.binning:cx-kwargs['y_gap_size_in_pixel']/2/self.binning+kwargs['y_gap_size_in_pixel']/self.binning] = 0
        if 'hole_diameter_in_pixel' in kwargs:
            if kwargs['hole_diameter_in_pixel'] != 0:
                cy = pylab.ceil((self.mask.shape[0]-1)/2.)
                cx = pylab.ceil((self.mask.shape[1]-1)/2.)
                X,Y = pylab.meshgrid(pylab.arange(0,self.mask.shape[1],1.0),
                                     pylab.arange(0,self.mask.shape[0],1.0))
                X -= cx
                Y -= cy
                R = pylab.sqrt(X**2 + Y**2)
                self.mask[R<=kwargs['hole_diameter_in_pixel']/(2.0*self.binning)] = 0
              
    def set_cy(self,cy=None):
        """
        Function sets vertical center position:
        =======================================
        Arguments:
        cy : vertical center position in pixel. If argument is None or not given center is set to the middle.
        """
        if cy == None:
            self.cy = (self.mask.shape[0]-1)/2.
        else:
            self.cy = cy

    def set_cx(self,cx=None):
        """
        Function sets horicontal center position:
        =========================================
        Arguments:
        cx : horizontal center position in pixel. If argument is None or not given center is set to the middle.
        """
        if cx == None:
            self.cx = (self.mask.shape[1]-1)/2.
        else:
            self.cx = cx

    def get_cx(self,option='unbinned'):
        if option == 'unbinned':
            return self.cx
        elif option == 'binned':
            if self.cx%1 == 0.5:
                return (self.cx-(self.binning-1)/2.)/(1.0*self.binning)
            else:
                return self.cx/(1.0*self.binning)
        else:
            print "ERROR: No valid option chosen."


    def get_cy(self,option='unbinned'):
        if option == 'unbinned':
            return self.cy
        elif option == 'binned':
            if self.cy%1 == 0.5:
                return (self.cy-(self.binning-1)/2.)/(1.0*self.binning)
            else:
                return self.cy/(1.0*self.binning)
        else:
            print "ERROR: No valid option chosen."
    
    def get_pixel_size(self,option='unbinned'):
        if option == 'unbinned':
            return self.pixel_size
        elif option == 'binned':
            return self.pixel_size*self.binning
        else:
            print "ERROR: No valid option chosen."

            
            
    # Functions that convert detector-coordinates:
    #def _get_q_from_r(self,r):
    #    return 4*pylab.pi/self._parent.source.photon.get_wavelength()*pylab.sin(pylab.arctan(r/self.distance)/2.0)
        #return 2.0*pylab.pi/self._parent.source.photon.get_wavelength()*r/self.distance
    
    #def _get_r_from_q(self,q):
    #    return self.distance*pylab.tan(pylab.arcsin(q*self._parent.source.photon.get_wavelength()/4.0/pylab.pi)*2.0)
        #return self.distance*q*self._parent.source.photon.get_wavelength()/2.0/pylab.pi
    
    #def _get_i_from_q(self,q):
    #    return int(round(self._get_r_from_q(self,q)/self.get_effective_pixelsize()))
   
    #def _get_q_from_i(self,i):
    #    return self._get_q_from_r(self,i*self.get_effective_pixelsize())        

