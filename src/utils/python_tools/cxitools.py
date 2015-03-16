#====================#
# Python tools - CXI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com


#from pylab import *
import h5py,os,numpy,time,sys

# CXI pixelmask bits
PIXEL_IS_PERFECT = 0
PIXEL_IS_INVALID = 1
PIXEL_IS_SATURATED = 2
PIXEL_IS_HOT = 4
PIXEL_IS_DEAD = 8
PIXEL_IS_SHADOWED = 16
PIXEL_IS_IN_PEAKMASK = 32
PIXEL_IS_TO_BE_IGNORED = 64
PIXEL_IS_BAD = 128
PIXEL_IS_OUT_OF_RESOLUTION_LIMITS = 256
PIXEL_IS_MISSING = 512
PIXEL_IS_IN_HALO = 1024
PIXEL_IS_ARTIFACT_CORRECTED = 2048
PIXEL_IS_IN_MASK = PIXEL_IS_INVALID |  PIXEL_IS_SATURATED | PIXEL_IS_HOT | PIXEL_IS_DEAD | PIXEL_IS_SHADOWED | PIXEL_IS_IN_PEAKMASK | PIXEL_IS_TO_BE_IGNORED | PIXEL_IS_BAD | PIXEL_IS_MISSING

CXIWRITER_DS_ACCUMULATIVE_PATTERN = -2
CXIWRITER_DS_SINGLE = -1

class CXIWriter:
    def __init__(self,filename,N,logger=None):
        self.filename = os.path.expandvars(filename)
        self.f = h5py.File(filename,"w")
        self.N = N
        self.logger = logger
    def write(self,d,prefix="",i=-1):
        for k in d.keys():
            name = prefix+"/"+k
            if self.logger is not None:
                self.logger.debug("Writing dataest %s",name)
            if isinstance(d[k],dict):
                if name not in self.f:
                    self.f.create_group(name)
                self.write(d[k],name,i)
            elif k != "i":
                self.write_to_dataset(name,d[k],d.get("i",i))
    def write_to_dataset(self,name,data,i):
        if self.logger is not None:
            self.logger.debug("Write dataset %s of event %i." % (name,i))
        if name not in self.f:
            #print name
            t0 = time.time()
            if numpy.isscalar(data):
                if i == -1:
                    s = [1]
                else:
                    s= [self.N]
                t=numpy.dtype(type(data))
                if t == "S":
                    t = h5py.new_vlen(str)
                axes = "experiment_identifier:value"
            else:
                data = numpy.array(data)
                s = list(data.shape)
                ndims = len(s)
                axes = "experiment_identifier"
                if ndims == 2: axes = axes + ":x"
                elif ndims == 3: axes = axes + ":y:x"
                elif ndims == 4: axes = axes + ":z:y:x"
                if i != -1:
                    s.insert(0,self.N)
                t=data.dtype
            self.f.create_dataset(name,s,t)
            self.f[name].attrs.modify("axes",[axes])
            t1 = time.time()
            if self.logger != None:
                self.logger.debug("Create dataset %s within %.1f sec.",name,t1-t0)

        if i == -1:
            if numpy.isscalar(data):
                self.f[name][0] = data
            else:
                self.f[name][:] = data[:]
        else:
            if numpy.isscalar(data):
                self.f[name][i] = data
            else:
                self.f[name][i,:] = data[:]
    def close(self):
        self.f.close()

class CXIReader:
    # location can be either a file or a directory
    def __init__(self,location,dsnames_stack={},**kwargs):
        self.logger = kwargs.get("logger",None)
        nevents = kwargs.get("nevents",0)
        ifirst = kwargs.get("ifirst",0)
        def_stack_ds = kwargs.get("def_stack_ds","/i")
        pick_events = kwargs.get("pick_events","in_sequence")
        event_filters = kwargs.get("event_filters",{})
        sorting_dsname = kwargs.get("sorting_dataset",None)
        self.dsnames_single = kwargs.get("dsnames_single",None)
        self.dsnames_stack = dsnames_stack

        [self.directories,self.filenames] = self._resolve_location(location)
        self.ifile = 0
        self.ifile_opened = None
        self.Nfiles = len(self.filenames)

        self.def_stack_ds = def_stack_ds
        self.Nevents_files = self.get_Nevents_files()
        self.Nevents_tot = 0
        for N in self.Nevents_files: self.Nevents_tot += N

        if self.logger != None:
            for d,f,N in zip(self.directories,self.filenames,self.Nevents_files):
                self.logger.info("Found file %s/%s with %i events.",d,f,N)

        self.ievent_file = -1
        self.ievent_tot = -1

        if nevents < 0:
            sys.exit("ERROR: Events to read smaller 0. Change your configuration.")
        elif nevents+ifirst > self.Nevents_tot:
            sys.exit("ERROR: Not enough events to read. Change your configuration.")
        
        to_process = []
        for N in self.Nevents_files: to_process.append(numpy.ones(N,dtype="bool"))
        to_process = self._filter_events(to_process,event_filters)
        to_process = self._pick_events(to_process,pick_events,ifirst,nevents)
        self.is_event_to_process = to_process
        self.ievent_process = -1
        self.Nevents_process = 0
        for N in to_process: self.Nevents_process += N.sum()

        self.order = self._set_order(to_process,sorting_dsname)

    def get_next(self):
        if self._next():
            return self._get(self.dsnames_stack,self.dsnames_single)
        else:
            return None

    def get_Nevents_files(self):
        N = []
        for i in range(len(self.filenames)):
            F = h5py.File(self.directories[i]+'/'+self.filenames[i],'r')
            N.append(F[self.def_stack_ds].shape[0])
            F.close()
        return N
        
    def close(self):
        if self.ifile_opened != None:
            self.F.close()

    def _resolve_location(self,location0):
        location = os.path.expandvars(location0)
        if os.path.isdir(location):
            fs = filter(lambda x: x[-4:] == ".cxi",os.listdir(location))
            fs.sort()
            directories = []
            filenames = []
            for f in fs:
                directories.append(location+"/")
                filenames.append(f.split("/")[-1])
        else:
            filenames = [location.split("/")[-1]]
            directories = [location[:-len(filenames[0])]]
        return [directories,filenames]        

    def _pick_events(self,to_process,mode,ifirst,nevents):
        temp = numpy.ones(self.Nevents_tot,dtype='bool')
        offset = 0
        for t in to_process:
            temp[offset:offset+len(t)] = t[:]        
            offset += len(t)
        if mode == 'in_sequence':
            temp[:ifirst] = False
            if nevents != 0:
                if nevents < temp.sum():
                    s = 0
                    for i in range(ifirst,self.Nevents_tot):
                        if temp[i]:
                            s += 1
                        if s == nevents:
                            break
                    temp[i+1:] = False
        elif mode == 'random':
            to_pick_from = arange(self.Nevents_tot)
            to_pick_from = list(to_pick_from[temp])
            temp = numpy.zeros_like(temp)
            for i in range(nevents):
                N = len(to_pick_from)
                if N != 1:
                    j = to_pick_from.pop(randint(N))
                    temp[j] = True
                else:
                    temp[to_pick_from[0]] = True
        else:
            print "ERROR: No valid picking mode chosen: %s" % mode
            return
        to_process_new = []
        i = 0
        for N in self.Nevents_files:
            to_process_new.append(temp[i:i+N])
            i += N
        return to_process_new

    def _filter_events(self,to_process,event_filters):
        for (i,dty,fle) in zip(range(self.Nfiles),self.directories,self.filenames):
            if self.logger != None:
                self.logger.info("Filtering %s/%s",dty,fle)
            f = h5py.File(dty+"/"+fle,"r")
            for flt_name in event_filters.keys():
                flt = event_filters[flt_name]
                filter_ds = f[flt["dataset_name"]].value.flatten()
                if "vmin" in flt.keys() and "vmax" in flt.keys():
                    F = (filter_ds >= flt["vmin"]) * (filter_ds <= flt["vmax"])
                elif "indices" in flt.keys():
                    if i != 0:
                        if self.logger != None:
                            self.logger.warning("Filter indices are applied to every file!")
                    F = numpy.zeros_like(to_process[i])
                    if not isinstance(flt["indices"],list):
                        F[filter_ds == flt["indices"]] = True
                    else:
                        for index in flt["indices"]:
                            F[filter_ds == index] = True
                else:
                    if self.logger != None:
                        self.logger.warning("No valid filter arguments given for filter %s!" % flt_name)
                    F = numpy.ones_like(to_process[i])
                to_process[i] *= F
                if self.logger != None:
                    self.logger.info("Filter %s - yield %.3f %% -> total yield %.3f %%",flt_name,100.*F.sum()/len(F),100.*to_process[i].sum()/len(F))
                    self.logger.info("Filter %s - First index: %i",flt_name,(arange(len(to_process[i]))[to_process[i]])[0])
        return to_process

    def _set_order(self,to_process,sorting_dsname):
        order = []
        for (i,Nevents_fle,dty,fle) in zip(range(self.Nfiles),self.Nevents_files,self.directories,self.filenames):
            if sorting_dsname == None:
                order.append(numpy.arange(Nevents_fle))
            else:
                order.append(numpy.zeros(Nevents_fle))
                if self.logger != None:
                    self.logger.debug("Set order of events in %s/%s",dty,fle)
                f = h5py.File(dty+"/"+fle,"r")
                sortdata = f[sorting_dsname]
                order[i][:] = argsort(sortdata)[:]
                f.close()
        return order


    # move to next event that shall be processed
    def _next(self):
        # skip events that shall not be processed
        while True:
            if self.ievent_process == self.Nevents_process:
                if self.logger != None:
                    self.logger.debug("Reached last event to process.")
                return False
            self.ievent_file += 1
            # return none if end of file list is reached
            if self.ifile >= self.Nfiles:
                if self.logger != None:
                    self.logger.debug("Reached end of list of files.")
                self.F.close()
                return False
            if self.ievent_file >= self.Nevents_files[self.ifile]:
                if self.logger != None:
                    self.logger.debug("Reached end of file (%i) %s/%s.",self.ifile,self.directories[self.ifile],self.filenames[self.ifile])
                self.ifile += 1
                self.ievent_file = -1
            if self.is_event_to_process[self.ifile][self.order[self.ifile][self.ievent_file]] == False:
                pass
                #if self.logger != None:
                #    self.logger.info("Skip event %i in file %i.",self.ievent_file,self.ifile)
            else:
                self.ievent_process += 1
                break
        return True

    def _open_file(self):
        if self.ifile_opened != self.ifile:
            if self.ifile_opened != None:
                if self.logger != None:
                    self.logger.debug("Closing file: %s/%s",self.directories[self.ifile_opened],self.filenames[self.ifile_opened])            
                self.F.close()
            if self.logger != None:
                self.logger.debug("Opening file: %s/%s",self.directories[self.ifile],self.filenames[self.ifile])
            self.F = h5py.File(self.directories[self.ifile]+'/'+self.filenames[self.ifile],'r')
            self.ifile_opened = self.ifile
        
    def _get(self,dsnames_stack,dsnames_single):
        self._open_file()
        D = {}
        D["i"] = self.ievent_process
        D["filename"] = self.filenames[self.ifile_opened]
        D["filename_i"] = self.ifile_opened
        if dsnames_stack != None:
            for (key,dsname) in dsnames_stack.items():
                self.logger.debug(dsname)
                D[key] = self.F[dsname][self.order[self.ifile][self.ievent_file]].copy()
        if dsnames_single != None:
            for (key,dsname) in dsnames_single.items():
                self.logger.debug(dsname)
                D[key] = self.F[dsname][:].copy()
        return D


def get_filters(C):
    filters = filter(lambda x: "filter" in x,C.keys())
    Cfilters = {}
    for f in filters:
        Cfilters[f] = C[f]
    return Cfilters


def intensities_cxi_to_h5(filename_input,filename_output,slice_i,ds_img="/imgXxX",ds_msk="/mskXxX",ds_cx=None,ds_cy=None,cropLength=None):
    f = h5py.File(filename_input,"r")
    cx = None
    cy = None
    if ds_cx != None:
        cx = f[ds_cx][slice_i]
    if ds_cy != None:
        cy = f[ds_cy][slice_i]
    img = f[ds_img][slice_i,:,:]
    img[img<0] = 0
    msk = f[ds_msk][slice_i,:,:]
    intensities_array_to_h5(img,msk,cx,cy,cropLength,filename_output)
    f.close()


# expects shifted pattern (center roughly in the middle)
def intensities_array_to_h5(img0,msk0,cx=None,cy=None,cropLength=None,save_to_file=None):
    import spimage,imgtools
    Nx = img0.shape[1]
    Ny = img0.shape[0]
    img = img0
    msk = msk0
    if cropLength != None:
        img = imgtools.crop(img0,cropLength,center=[cy,cx])
        msk = imgtools.crop(msk0,cropLength,center=[cy,cx])
        Nx = cropLength
        Ny = cropLength
    elif cx != None and cy != None:
        img = imgtools.recenter(img0,cx,cy)
        msk = imgtools.recenter(msk0,cx,cy)
    I = spimage.sp_image_alloc(Nx,Ny,1)
    I.image[:,:] = img[:,:]
    I.mask[:,:] = msk[:,:]
    I2 = spimage.sp_image_shift(I)
    I2.shifted = 0
    I2.phased = 0
    spimage.sp_image_free(I)
    if save_to_file != None:
        spimage.sp_image_write(I2,save_to_file,0)
    else:
        return I2


def write_sigma_map(filename_cheetah,scale=1.):
    fc = h5py.File(filename_cheetah,"r+")
    M = fc["/entry_1/image_2/mask"]
    H = fc["/cheetah/shared/detector0_class0_sigma_downsampled"]
    electronics_sigma = 2.5
    artifact_sigma = 50.
    N = min([M.shape[1],M.shape[2]])
    try:
        fc.create_group("/sigma_map")
    except:
        print "Group already exists"
    s = (M.shape[0],N,N)
    t = "float"
    axes = "experiment_identifier:y:x"
    try:
        fc.create_dataset("/sigma_map/data",s,t)
        fc["/sigma_map/data"].attrs.modify("axes",axes)
    except:
        print "Dataset already exists"
    for i in arange(s[0]):
        fc["/sigma_map/data"][i,:,:] = (electronics_sigma+H[:N,:N]+((M[i,:N,:N] & PIXEL_IS_ARTIFACT_CORRECTED) != 0)*artifact_sigma)*scale
    fc.close()




