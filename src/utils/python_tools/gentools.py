#========================#
# Python tools - general #
#========================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

import numpy,os,re,csv,string,time,datetime,logging
logger = logging.getLogger("pythontools")

# in_filter can be a string or a list of strings
def get_filenames(in_filter=None,path="./"):
    filenames = os.popen('ls %s' % path).readlines()
    if in_filter:
        if isinstance(in_filter,list):
            for filt in in_filter:
                filenames = [f for f in filenames if re.search(filt,f)]
        else:
            filenames = [f for f in filenames if re.search(in_filter,f)]
    filenames.sort()
    return filenames

def estimate_type(var):
    #first test bools
    if var.lower() == 'true':
        return True
    elif var.lower() == 'false':
        return False
    elif var.lower() == 'none':
        return None
    else:
        #int
        try:
            return int(var)
        except ValueError:
            pass
        #float
        try:
            return float(var)
        except ValueError:
            pass
        #string
        try:
            return str(var)
        except ValueError:
            raise NameError('Something messed up autocasting var %s (%s)' % (var, type(var)))

#def png_mask_to_h5(filename,Nx,Ny):
#    import Image,spimage
#    I = Image.open(filename)
#    D = numpy.array(I.getdata())[:,0]
#    D=D.reshape((Nx,Ny))
#    img = spimage.sp_image_alloc(Nx,Ny,1)
#    img.mask[:,:] = D[:,:]/255
#    spimage.sp_image_write(img,filename[:-4]+'.h5',0)
#    spimage.sp_image_free(img)

def get_png_mask(filename):
    try: 
        from pillow import Image
    except:
        try:
            import Image
        except:
            import image as Image
    I = Image.open(filename)
    D = numpy.array(I.getdata())
    if len(D.shape) > 1:
        D = numpy.array(I.getdata())[:,0]
    D = D.reshape((I.size[1],I.size[0]))
    D = D[:,:]/255.
    D = D.round()
    D = numpy.int16(D)
    return D


def get_png_grayvalue(filename):
    try: 
        from pillow import Image
    except:
        import Image
    I = Image.open(filename)
    D = numpy.array(I.getdata())
    if len(D.shape) > 1:
        D = numpy.array(I.getdata())[:,0]
    D = D.reshape((I.size[1],I.size[0]))
    return D

def save_to_csv(filename,list_of_arrays,list_of_array_names=[]):
    """Save given array values to a new csv file."""
    f = open(filename, 'wb')
    writer = csv.writer(f)
    if list_of_array_names != []: writer.writerow(list_of_array_names)
    for i in range(0,len(list_of_arrays[0])):
        row = []
        for j in range(0,len(list_of_arrays)): 
            row.append(list_of_arrays[j][i])
        writer.writerow(row)
    f.close()

def load_from_csv(filename,has_header=False):
    """Loads csv file and gives out one array for each column and the header."""
    f = open(filename,'r')
    reader = csv.reader(f)
    rows = []
    for row in reader: rows.append(row)
    f.close()
    if has_header:
        header = rows[0]
        rows.remove(rows[0])
    else:
        header = []
        for i in range(0,len(row)):
            header.append(string.letters[i])
    print header
    out = []
    for i in range(0,len(row)): out.append([header[i],[]])
    out = dict(out)
    for i in range(0,len(rows)):
        for j in range(0,len(row)):
            out[header[j]].append(estimate_type(rows[i][j]))
    return out

def create_numberstring(number,number_of_letters):
    number_str = ""
    for j in -numpy.arange(-number_of_letters+1,1,1):
        number_str =  number_str + ("%i" % (number/pow(10,j)%10))
    return number_str


def smoothList(list,strippedXs=False,degree=10):  
    if strippedXs==True: return Xs[0:-(len(list)-(len(list)-degree+1))]  
    smoothed=[0]*(len(list)-degree+1)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(list[i:i+degree])/float(degree)  
    return smoothed  

def smoothListTriangle(list,strippedXs=False,degree=5):  
    weight=[]  
    window=degree*2-1  
    smoothed=[0.0]*(len(list)-window)  
    for x in range(1,2*degree):weight.append(degree-abs(degree-x))  
    w=numpy.array(weight)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*w)/float(sum(w))  
    return smoothed  

def smoothListGaussian(list,strippedXs=False,degree=5):  
    window=degree*2-1  
    weight=numpy.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(numpy.exp((4*(frac))**2))  
        weightGauss.append(gauss)  
    weight=numpy.array(weightGauss)*weight  
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
    return smoothed 

def smoothGaussian(array,sigma):
    return numpy.ifft(numpy.fft(array)*1/numpy.sqrt(2*numpy.pi)/sigma*numpy.exp(numpy.arange(0,len(array),1.0)**2/2.0/sigma**2))

def smoothGaussian1dMirror(array,sigma):
    array_reversed = list(array.copy())
    array_reversed.reverse()
    array_reversed = numpy.array(array_reversed)
    array_mirrored = numpy.zeros(2*len(array))
    array_mirrored[:len(array)] = array[:]
    array_mirrored[len(array):] = array_reversed[:]
    array_smoothed = numpy.ifft(numpy.fft(array_mirrored)*1/numpy.sqrt(2*numpy.pi)/sigma*numpy.exp(numpy.arange(0,len(array_mirrored),1.0)**2/2.0/sigma**2) )
    array_smoothed = array_smoothed[:len(array)]
    print array_smoothed
    return array_smoothed

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
         
    input:
    x: the input signal 
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.
   
    output:
    the smoothed signal
           
    example:
    
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    
    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    
    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
       #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
        
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1] 

gaussian = lambda x,A0,x0,sigma: A0*numpy.exp((-(x-x0)**2)/(2.*sigma**2))
double_gaussian = lambda x,p1,p2: gaussian(x,p1[0],p1[1],p1[2])+gaussian(x,p2[0],p2[1],p2[2])

def double_gaussian_fit(xdata=None,ydata=None,p_init=None,show=False):
    from scipy.optimize import leastsq

    if xdata == None and ydata == None:
        # generate test data
        A01,A02 = 7.5,10.4
        mean1, mean2 = -5, 4.
        std1, std2 = 7.5, 4. 
        xdata =  numpy.linspace(-20, 20, 500)
        ydata = double_gaussian(xdata,[A01,mean1,std1],[A02,mean2,std2])

    if p_init == None:
        p_A01, p_A02, p_mean1, p_mean2, p_sd1, p_sd2  = [ydata.max(),
                                                         ydata.max(),
                                                         xdata[len(xdata)/3],
                                                         xdata[2*len(xdata)/3],
                                                         (xdata.max()-xdata.min())/10.,
                                                         (xdata.max()-xdata.min())/10.]
        p_init = [p_A01, p_mean1, p_sd1,p_A02, p_mean2, p_sd2] # Initial guesses for leastsq
    else:
        [p_A01, p_mean1, p_sd1,p_A02, p_mean2, p_sd2] = p_init # Initial guesses for leastsq

    err = lambda p,x,y: abs(y-double_gaussian(x,[p[0],p[1],p[2]],[p[3],p[4],p[5]]))

    plsq = leastsq(err, p_init, args = (xdata, ydata))

    p_result = [[A1, mean1, s1], [A2, mean2, s2]] = [[plsq[0][0],plsq[0][1],plsq[0][2]],[plsq[0][3],plsq[0][4],plsq[0][5]]]
    yest = double_gaussian(xdata,p_result[0],p_result[1])

    if show:
        import pylab
        pylab.figure()
        yinit = double_gaussian(xdata,[p_A01,p_mean1,p_sd1],[p_A02,p_mean2,p_sd2])
        yest1 = gaussian(xdata,plsq[0][0],plsq[0][1],plsq[0][2])
        yest2 = gaussian(xdata,plsq[0][3],plsq[0][4],plsq[0][5])
        pylab.plot(xdata, ydata, 'r.',color='red', label='Data')
        pylab.plot(xdata, yinit, 'r.',color='blue', label='Starting Guess')
        pylab.plot(xdata, yest, '-',lw=3.,color='black', label='Fitted curve (2 gaussians)')
        pylab.plot(xdata, yest1, '--',lw=1.,color='black', label='Fitted curve (2 gaussians)')
        pylab.plot(xdata, yest2, '--',lw=1.,color='black', label='Fitted curve (2 gaussians)')
        pylab.legend()
        pylab.show()
        
    return [p_result, yest]

def gaussian_fit(xdata=None,ydata=None,p_init=None,show=False):
    from scipy.optimize import leastsq
    
    if xdata == None and ydata == None:
        # Generate test data
        A0_test = 5.5
        x0_test = -2.
        s_test = 1. 
        xdata =  numpy.linspace(-20., 20., 500)
        ydata = gaussian(xdata,A0_test,x0_test,s_test)

    if p_init == None:
        # Guess parameters
        A0_init, x0_init, s_init  = [ydata.max(),
                                     xdata[len(xdata)/2].mean(),
                                     (xdata.max()-xdata.min())/10.]
        p_init = [A0_init, x0_init, s_init]
    else:
        [A0_init, x0_init, s_init] = p_init

    err = lambda p,x,y: abs(y-gaussian(x,p[0],p[1],p[2]))
    plsq = leastsq(err, p_init, args = (xdata, ydata))

    p_result = [A0_result, x0_result, s_result] = [plsq[0][0],plsq[0][1],plsq[0][2]]
    yest = gaussian(xdata,A0_result, x0_result, s_result)

    if show:
        import pylab
        pylab.figure()
        yinit = gaussian(xdata,A0_init, x0_init, s_init)
        pylab.plot(xdata, ydata, 'r.',color='red', label='Data')
        pylab.plot(xdata, yinit, 'r.',color='blue', label='Starting Guess')
        pylab.plot(xdata, yest, '-',lw=3.,color='black', label='Fitted curve')
        pylab.legend()
        pylab.show()
        
    return [p_result, yest]

def bootstrap_gaussian_fit(xdata,ydata,p_init0=None,show=False,Nfract=0.5,n=100):
    if p_init0 == None:
        p_init = gaussian_fit(xdata,ydata)[0]
    else:
        p_init = p_init0
    ps = []
    N = round(len(xdata)*Nfract)
    for i in range(n):
        random_pick = numpy.random.randint(0,N,(N,))
        xdata1 = xdata[random_pick]
        ydata1 = ydata[random_pick]
        variation = 0.5
        p0 = numpy.array(p_init) * (1+((-0.5+numpy.random.rand(len(p_init)))*variation))
        ps.append(gaussian_fit(xdata1,ydata1,tuple(p0),show)[0])
    ps = numpy.array(ps)
    p_result = ps.mean(0)
    p_std = ps.std(0)
    yest = gaussian(xdata,p_result[0], p_result[1], p_result[2])
    return [p_result, yest, p_std]
            


def my_imsave(fname, arr, **kwargs): imsave(fname, arr, **kwargs)

def imsave(fname, arr, **kwargs):

    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    from matplotlib.colors import ColorConverter as CC
    C = CC()
    
    import pylab
    if pylab.isinteractive():
        i_was_on = True
        pylab.ioff()
    else:
        i_was_on = False
        
    fig = Figure(figsize=arr.shape[::-1], dpi=1, frameon=False)
    canvas = FigureCanvas(fig)

    if 'background' in kwargs.keys():
        if kwargs['background'] != 'transparent':
            [r,g,b] = C.to_rgb(kwargs['background'])
            BG = numpy.ones(shape=(arr.shape[0],arr.shape[1],3))
            BG[:,:,0] = r
            BG[:,:,1] = g
            BG[:,:,2] = b
            fig.figimage(BG)

    fig.figimage(arr,
                 xo = kwargs.get('xo',0),
                 yo = kwargs.get('yo',0),
                 alpha = kwargs.get('alpha',None),
                 norm = kwargs.get('norm',None),
                 cmap = kwargs.get('cmap',None),
                 vmin = kwargs.get('vmin',None),
                 vmax = kwargs.get('vmax',None),
                 origin = kwargs.get('origin',None))
    
    fig.savefig(fname,
                dpi=1,
                format = kwargs.get('format',None))
    if i_was_on:
        pylab.ion()

def dict_to_dict(dict_src,dict_dest):
    for k in dict_src.keys():
        dict_dest[k] = dict_src[k]
        
def read_configfile(configfile):
    import ConfigParser
    config = ConfigParser.ConfigParser()
    with open(configfile,"r") as f:
        config.readfp(f)
        confDict = {}
        for section in config.sections(): 
            confDict[section] = {}
            c = config.items(section)
            for (key,value) in c:
                v = estimate_type(value)
                if isinstance(v,str):
                    if "[" == v[0] and "]" == v[-1]:
                        v = v[1:-1].split(",")
                        v = [w for w in v if len(w) > 0]
                        for i in range(len(v)):
                            if '$' in v[i]:
                                v[i] = os.path.expandvars(v[i])
                            v[i] = estimate_type(v[i]) 
                    else:
                        if '$' in v:
                            v = os.path.expandvars(v)
                confDict[section][key] = v
    return confDict

class Configuration:
    def __init__(self,config={},default={},verbose=False):
        self.verbose = verbose
        if isinstance(config,str):
            self.confDict = read_configfile(config)
            self.configfile = config
        else:
            self.confDict = config
        
        if isinstance(default,str):
            defDict = read_configfile(default)
        else:
            defDict = default
        self.set_unspecified_to_default(defDict)
    
    def set_unspecified_to_default(self,defaultDict):
        for section in defaultDict.keys():
            matched_sections = [s for s in self.confDict.keys() if section in s]
            if matched_sections == []:
                self.confDict[section] = {}
                logger.info("Add section %s to configuration as it did not exist." % section)
            for variableName in defaultDict[section].keys():
                for ms in matched_sections:
                    if variableName not in self.confDict[ms].keys():
                        self.confDict[ms][variableName] = defaultDict[section][variableName]
                        logger.info("Add variable %s with default value %s to configuration section %s as variable did not exist." % (variableName,str(defaultDict[section][variableName]),ms))

    def write_to_file(self,filename):
        ls = ["# Configuration file\n# Automatically written by Configuration instance\n\n"]
        for section_name,section in self.confDict.items():
            if isinstance(section,dict):
                ls.append("[%s]\n" % section_name)
                for variable_name,variable in section.items():
                    ls.append("%s=%s\n" % (variable_name,str(variable)))
                ls.append("\n")
        s = open(filename,"w")
        s.writelines(ls)
        s.close()        

def mkdir_timestamped(dirname):
    ts = time.time()
    time_string = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M%S')
    real_dir = dirname+'_'+time_string
    os.system('mkdir %s' % real_dir)
    return real_dir

def make_link(linkname,targetname):
    if linkname[-1] == "/":
        ln = linkname[:-1]
    else:
        ln = linkname
    os.system('rm -f %s' % linkname)
    os.system('ln -s %s %s' % (targetname,linkname))
