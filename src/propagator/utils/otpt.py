

class Output:
    """
    OUTPUT of propagator provides user with results and functions for plotting.
    """
    def __init__(self,amplitudes,input_object):
        self.amplitudes = amplitudes.copy()
        self.input_object = input_object 
    
    def get_intensity_pattern(self,mask_cut=False):
        """
        Returns 2-dimensional array with intensity values in photons per pixel (binned).
        """
        I = abs(self.amplitudes)**2
        if mask_cut:
            [X_min,X_max,Y_min,Y_max] = self._get_pattern_limits()
            return I[Y_min:Y_max+1,X_min:X_max+1]
        else:
            return I

    def get_intensity_radial_average(self):
        """
        Returns 1-dimensional array with intensity average in photons per pixel (binned). x-coordinate sampling is pixel (binned). 
        """
        I = self.get_intensity_pattern()
        return imgutils.radial_pixel_average(I,
                                             cx=self.input_object.detector.get_cx('binned'),
                                             cy=self.input_object.detector.get_cy('binned'))


    def get_intensity_radial_sum(self):
        """
        Returns 1-dimensional array with intensity average in photons per pixel (binned). x-coordinate sampling is pixel (binned). 
        """
        I = self.get_intensity_pattern()
        return imgutils.radial_pixel_sum(I,
                                         cx=self.input_object.detector.get_cx('binned'),
                                         cy=self.input_object.detector.get_cy('binned'))

            
    def plot_radial_distribution(self,scaling="binned pixel and nyquist pixel",mode="all",noise=None):
        """
        Creates 1-dimensional plot(s) showing radial distribution of scattered photons.
        Usage: plot_radial_distribution([scaling],[mode],[noise])
        Arguments:
        - scaling: Specifies spatial scaling.
                   Can be set to 'binned pixel', 'nyquist pixel', 'binned pixel and nyquist pixel' or 'meter'.
                   'binned pixel and nyquist pixel' leads to creation of two plots in one figure using pixel- and Nyquist-pixel-scaling.
        - mode:    Mode specifies whether the radial average or the radial sum will be plotted.
                   Can be set to 'radial average', 'radial sum' or 'all'.
        - noise:   Specifies noise and can be set to 'poisson'.
        """
        Ir_avg = self.get_intensity_radial_average()
        Ir_sum = self.get_intensity_radial_sum()
        if noise == 'poisson':
            def noise(data): return pylab.poisson(data)
        else:
            def noise(data): return data
        def get_arguments(sc):
            if mode == "all":
                legend_args = [('Radial sum', 'Radial average'),'upper right']
                if sc == "binned pixel":
                    r = numpy.arange(0,len(Ir_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
                plot_args = [r,noise(Ir_sum),'k',r,noise(Ir_avg),'k:']
            else:
                if sc == "binned pixel":
                    r = numpy.arange(0,len(Ir_sum),1)
                elif sc == "nyquist pixel":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
                elif sc == "meter":
                    r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2*self.pixel_size/len(Ir_sum))
                if mode == "radial sum":
                    legend_args = [('Radial sum'),'upper right']
                    plot_args = [r,noise(Ir_sum),'k']
                elif mode == "radial average":
                    legend_args = [('Radial average'),'upper right']
                    plot_args = [r,noise(Ir_avg),'k']
            return [plot_args,legend_args]

        if scaling == "binned pixel and nyquist pixel":
            f1d = pylab.figure(figsize=(10,5))
            f1d.suptitle("\nRadial distribution of scattered photons in detector plane", fontsize=16)
            str_scaling = "binned pixel"
            f1d_ax_left = f1d.add_axes([0.1, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            str_scaling = "nyquist pixel"
            f1d_ax_right = f1d.add_axes([0.55, 0.1, 0.35, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
            [plot_args,legend_args] = get_arguments('binned pixel')
            f1d_ax_left.semilogy(*plot_args)
            f1d_ax_left.legend(*legend_args)
            [plot_args,legend_args] = get_arguments('nyquist pixel')
            f1d_ax_right.semilogy(*plot_args)
            f1d_ax_right.legend(*legend_args)
            f1d.show()
            return
        elif scaling == "binned pixel":
            str_scaling = "binned pixel"
            r = numpy.arange(0,len(Ir_sum),1)
        elif scaling == "nyquist pixel":
            str_scaling == "nyquist pixel"
            r = numpy.arange(0,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2,min([self.nyquistpixel_number_x,self.nyquistpixel_number_y])/2/len(Ir_sum))
        elif scaling == "meter":
            str_scaling = "meter"
            r = numpy.arange(0,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size,min([self.pixel_number_x,self.pixel_number_y])/2*self.pixel_size/len(Ir_sum))
        else:
            print "ERROR: %s is no valid scaling" % scaling
            return
        [plot_args,legend_args] = get_arguments(r,scaling)
        f1d = pylab.figure(figsize=(5,5))
        f1d.suptitle("\nRadial distribution of scattered photons in detector plane", fontsize=16)
        f1d_ax = f1d.add_axes([0.2, 0.1, 0.7, 0.7],title='Radial scaling:' + str_scaling,xlabel="r [" + str_scaling + "]",ylabel="I(r) [photons/" + str_scaling + "]")
        f1d_ax.semilogy(*plot_args)
        f1d_ax.legend(*legend_args)
        f1d.show()

    def get_nyquist_pixel_size(self):
        return proptools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),self.input_object.sample.get_area())

    def _get_gapsize(self,X_min,X_max,Y_min,Y_max):
        """
        Returns gapsize of pattern in pixels (binned)
        """
        gapsize = 0
        M = self.input_object.detector.mask
        for i in pylab.arange(X_min,X_max+1,1):
            if (M[:,i]==0).all():
                for j in pylab.arange(i,X_max+1,1):
                    if (M[:,j]==1).any():
                        gapsize = j-i
                        break
                break
        for i in pylab.arange(Y_min,Y_max+1,1):
            if (M[i,:]==0).all():
                for j in pylab.arange(i,Y_max+1,1):
                    if (M[j,:]==1).any():
                        gapsize = j-i
                        break
                break
        return gapsize

    def _get_pattern_limits(self):
        """
        Returns spatial limits of pattern in pixels (binned)
        """
        X_min = 0
        Y_min = 0
        X_max = self.amplitudes.shape[1]
        Y_max = self.amplitudes.shape[0]
        M = self.input_object.detector.mask
        for i in pylab.arange(0,M.shape[1],1):
            if (M[:,i]==1).any():
                X_min = i
                break
        for i in M.shape[1]-pylab.arange(1,M.shape[1],1):
            if (M[:,i]==1).any():
                X_max = i
                break
        for i in pylab.arange(0,M.shape[0],1):
            if (M[i,:]==1).any():
                Y_min = i
                break
        for i in M.shape[0]-pylab.arange(1,M.shape[0],1):
            if (M[i,:]==1).any():
                Y_max = i
                break
        return [X_min,X_max,Y_min,Y_max]

    def plot_pattern(self,**kwargs):
        """
        Function plots the 2d intensity pattern of scattered photons.
        =============================================================

        Keyword arguments:

        - scaling:          'nyquist', 'meter', 'pixel' or 'binned pixel' (default)
        - noise:            'poisson' or 'none' (default)
        - logscale:         False / True (default)
        - saturation_level:  False (default) / True
        - use_gapmask:      False / True (default)
        - dpi:              dpi (default 200)

        """

        from matplotlib.colors import LogNorm

        if self.input_object.detector.binning == 1:
            scaling = 'pixel'
        else:
            scaling = 'binned pixel'
        scalingargs = ['nyquist','meter','binned pixel','pixel']
        noise = 'none'
        noiseargs = ['poisson','none']
        logscale = True
        logscaleargs = [False,True]
        saturationlevel = False
        saturationlevelargs = [False,True]
        use_gapmask = True
        use_gapmaskargs = [False,True]

        if 'dpi' not in kwargs.keys():
            kwargs['dpi'] = 200
        if 'outfile' not in kwargs.keys():
            kwargs['outfile'] = False
        if 'bg' not in kwargs.keys():
            kwargs['bg'] = 0.0
        I = self.get_intensity_pattern()

        optionkeys = ["scaling","noise","logscale","saturation_level","use_gapmask","outfile","dpi","bg"]
        options = [scaling,noise,logscale,saturationlevel,use_gapmask]
        optionargs = [scalingargs,noiseargs,logscaleargs,saturationlevelargs,use_gapmaskargs,[],[],[]]
        keys = kwargs.keys()
        for i in range(0,len(keys)):
            key = keys[i]
            if not key in optionkeys:
                print "ERROR: %s is not a proper key." % key
                return
            keyarg = kwargs[key]
            j = optionkeys.index(key)
            if not keyarg in optionargs[j]:
                err = False
                if key == 'dpi':
                    if not pylab.isreal(keyarg):
                        err = True
                elif key == 'outfile':
                    if keyarg != False:
                        if keyarg[-4:] != '.png':
                            err = True
                elif key == 'bg':
                    if keyarg>1. or keyarg<0.:
                        err = True
                if err:
                    print "ERROR: %s is not a proper argument for %s." % (keyarg,key)
                    return
            exec "%s = '%s'" % (key,keyarg)
        
        eff_pixel_size_detector = self.input_object.detector.pixel_size * self.input_object.detector.binning
        pixel_size_detector = self.input_object.detector.pixel_size
        pixel_size_nyquist = proptools.get_nyquist_pixel_size(self.input_object.detector.distance,self.input_object.source.photon.get_wavelength(),self.input_object.sample.get_area())
        if scaling == "nyquist":
            I *= eff_pixel_size_nyquist**2/eff_pixel_size_detector**2
            u = eff_pixel_size_detector/pixel_size_nyquist
            str_scaling = "Nyquist pixel"
        elif scaling == "meter":
            I /= eff_pixel_size_detector**2
            u = eff_pixel_size_detector
            str_scaling = "m^2"
        elif scaling == "pixel":
            u = self.input_object.detector.binning
            str_scaling = scaling
        elif scaling == "binned pixel":
            u = 1.0
            str_scaling = scaling

        I /= u**2

        if noise == "poisson":
            I = abs(pylab.poisson(I))

        if saturationlevel and self.input_object.detector.saturationlevel > 0:
            I *= u**2/(1.0*self.input_object.detector.binning**2)
            I[I>self.input_object.detector.saturationlevel] = self.input_object.detector.saturationlevel
            I /= u**2/(1.0*self.input_object.detector.binning**2)

        [X_min,X_max,Y_min,Y_max] = self._get_pattern_limits()
        xlimit = u*(X_max-X_min)
        ylimit = u*(Y_max-Y_min)
        gapsize = self._get_gapsize(X_min,X_max,Y_min,Y_max)

        I = I[Y_min:Y_max+1,X_min:X_max+1]

        if str_scaling == "binned pixel":
            if self.input_object.detector.binning == 1:
                str_scaling_label = "pixel"
            else:
                str_scaling_label = "%ix%i binned pixel" % (self.input_object.detector.binning,self.input_object.detector.binning)
        else:
            str_scaling_label = str_scaling

        if use_gapmask:
            M = self.input_object.detector.mask.copy()
            if saturationlevel and self.input_object.detector.saturationlevel > 0:
                M[I>=self.input_object.detector.saturationlevel] = 0
            M[M==0] *= pylab.nan

        if logscale:
            #I = pylab.log10(I)
            I *= M
            #str_Iscaling = r"$\log_{10}\left( I \left[ \frac{\mathrm{photons}}{\mathrm{%s}} \right] \right)$" % str_scaling
            str_Iscaling = r"Number of photons"
        else:
            str_Iscaling = r"$I \left[ \frac{\mathrm{photons}}{\mathrm{%s}} \right]$" % str_scaling

        Wsizey = 9
        Wsizex = 9
        fsize = 12
        pylab.clf()
        fig = mpy.figure(1,figsize=(Wsizex,Wsizey))
        mpy.draw()
        fig.set_size_inches(Wsizex,Wsizey)
        fig.suptitle("- PROPAGATOR -", fontsize=fsize+2)
        alignment = {'horizontalalignment':'center','verticalalignment':'center'}

        fig.text(0.5,(16.75/18.0),r"$E_{\mathrm{photon}} = %.0f$ eV ; $\lambda = %.2f$ nm ; $I = %.1e$ $\frac{\mathrm{ph}}{\mu\mathrm{m}^2}$ ; $N_{\mathrm{ph,sc}} = %.1e$ ph ; $D_{\mathrm{detector}} = %0.3f$ mm" %  (self.input_object.source.photon.get_energy("eV"),
                                                                                                                                                                    self.input_object.source.photon.get_wavelength()/1.0E-09,
                                                                                                                                                                    self.input_object.source.pulse_energy/self.input_object.source.photon.get_energy()/(self.input_object.source.get_area()/1E-12),
                                                                                                                                                                    self.get_intensity_pattern().sum(),
                                                                                                                                                                    self.input_object.detector.distance/1.0E-03),
                 fontsize=fsize,bbox=dict(fc='0.9',ec="0.9",linewidth=10.0),**alignment) 

        ax = fig.add_axes([3/15.0,5/18.0,10/15.0,10/18.0])
        ax.set_xlabel(r"$x$ [" + str_scaling_label + "]",fontsize=fsize)
        ax.xaxis.set_label_position('top')
        ax.set_ylabel(r"$y$ [" + str_scaling_label + "]",fontsize=fsize)

        axcolor = fig.add_axes([3/15.0,3.5/18.0,10/15.0,0.5/18.0])
        for a in [ax,axcolor]:
            for label in a.xaxis.get_ticklabels():
                label.set_fontsize(fsize)
            for label in a.yaxis.get_ticklabels():
                label.set_fontsize(fsize)

        if logscale:
            bgcol = kwargs['bg']
            imv = ax.matshow(bgcol*pylab.ones_like(I),extent=[-xlimit/2,xlimit/2,-ylimit/2,ylimit/2],interpolation="nearest",cmap='gray',vmax=1.,vmin=0.)
            im = ax.matshow(I,extent=[-xlimit/2,xlimit/2,-ylimit/2,ylimit/2],interpolation="nearest",norm=LogNorm())
        else:
            im = ax.matshow(I,extent=[-xlimit/2,xlimit/2,-ylimit/2,ylimit/2],interpolation="nearest",)
        cb = fig.colorbar(im, cax=axcolor,orientation='horizontal')
        cb.set_label(str_Iscaling,fontsize=fsize)

        oversampling_ratio = pixel_size_nyquist/eff_pixel_size_detector
        oversampling_ratio_wo_binning = pixel_size_nyquist/pixel_size_detector
        D = self.input_object.detector.distance
        A =  self.input_object.sample.get_area()
        wavelength = self.input_object.source.photon.get_wavelength()
        [res_horizontally,res_vertically] = self.input_object.get_max_achievable_crystallographic_resolution()
        res_corner = 1/pylab.sqrt(1/res_horizontally**2 + 1/res_vertically**2)
        miss_Ny = gapsize*eff_pixel_size_detector/pixel_size_nyquist
        fig.text(0.5,(1./18.0),r"\textbf{Properties}\\ Linear oversampling ratio: $%.2f$ (binning $%i\times%i$) ; $%.2f$ (no pixel binning)\\" % (oversampling_ratio,self.input_object.detector.binning,self.input_object.detector.binning,oversampling_ratio_wo_binning)+
                 r"Crystallographic resolution (full period): $%.1f$ nm (horizontal) ; $%.1f$ nm (vertical) ; $%.1f$ nm (corner)\\" % (res_horizontally/1.0E-09,res_vertically/1.0E-09,res_corner/1.0E-09)+
                 r"Gap width: $g=%.2f\mathrm{ mm}=%.1f$ Nyquist pixels" % (gapsize*eff_pixel_size_detector/1.0E-03,miss_Ny),fontsize=fsize,bbox=dict(fc='0.9',ec="0.9",linewidth=10.0),**alignment)
        #if miss_Ny>2.8:
        #    print "\n!!!\nMissing mode(s) expected (gap width: %.1f Nyquist pixels) \n\nTweaking of one of the parameters recommended:\n- Wavelength w = %.2f nm\n- Sample radius r = %.0f nm\n- Gap size g = %.1f mm\n- Detector distance d = %.0f mm" % (miss_Ny,(rec_wavelength+0.01E-9)*1.0E9,(rec_r-1.0E-9)*1.0E9,(rec_gapsize-0.1E-3)*1.0E3,(rec_d+1.0E-3)*1.0E3)

        if kwargs['outfile'] == False:
            pylab.show()
        else:
            mpy.savefig(kwargs['outfile'],dpi=kwargs['dpi'])
   
    def save_pattern_to_file(self,filename,**kwargs):
        """
        Function saves dataset to file of specified format.
        ===================================================
        
        Arguments:
        - filename: The file-format is specified using one of the following file-endings:
                    - \'.h5\'
                    - \'.png\'

        
        Keyword arguments:
        - log: True / False (default)
        - poisson: True / False (default)
        - colorscale: \'jet\' (default) / \'gray\'
        - use_spimage: True / False (default)

        """
        pattern = self.get_intensity_pattern()
        mask = self.input_object.detector.mask
        if 'poisson' in kwargs:
            if kwargs['poisson']:
                pattern = pylab.poisson(pattern)
        if 'log' in kwargs:
            if kwargs['log']:
                pattern = pylab.log10(pattern)
                pattern[pylab.isfinite(pattern)==False] = pattern[pylab.isfinite(pattern)].min()
        use_spimage = kwargs.get('use_spimage',False)
        colorscale = kwargs.get('colorscale','jet')
        if use_spimage:
            import spimage
            if filename[-3:]=='.h5':
                color = 0
            elif filename[-3:]=='.png':
                if colorscale  == 'gray':
                        color = 1
                elif colorscale  == 'jet':
                        color = 16
            else:
                print "ERROR: %s is not a valid fileformat for this function." % filename[-3:]
                return
            tmp_data = spimage.sp_image_alloc(pattern.shape[1],pattern.shape[0],1)
            tmp_data.image[:,:] = pattern[:,:]
            tmp_data.mask[:,:] = mask[:,:]
            spimage.sp_image_write(tmp_data,filename,0)
            spimage.sp_image_free(tmp_data)
        else:
            if filename[-4:]=='.png':
                pylab.imsave(filename,pattern*pylab.log10(mask*10),cmap=pylab.cm.get_cmap(colorscale))
            elif filename[-3:]=='.h5':
                import h5py
                f = h5py.File(filename,'w')
                pattern_ds = f.create_dataset('intensities', pattern.shape, pattern.dtype)
                pattern_ds[:,:] = pattern[:,:]
                amplitudes_ds = f.create_dataset('amplitudes', self.amplitudes.shape, self.amplitudes.dtype)
                amplitudes_ds[:,:] = self.amplitudes[:,:]
                f.close()

