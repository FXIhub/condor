


class Material:
    """
    A class of a sample object.
    Sample material.

    """
    def __init__(self,parent=None,**kwargs):
        # Check for valid keyword arguments
        allk = ["massdensity","material_type"]
        self._unproc_kws = [k for k in kwargs.keys() if (k not in allk and k[0] != "c" and len(k) > 3)]
        if len(self._unproc_kws) > 0:
            print self._unproc_kws
            logger.error("Material object initialisation failed due to illegal keyword arguments.")
            return
        # Start intialisation
        self._parent = parent
        if "massdensity" in kwargs:
            self.materialtype = "custom"
            for key in kwargs:
                if key[0] == 'c' or key == 'massdensity':
                    exec "self." + key + " = kwargs[key]"
                else:
                    logger.error("%s is no valid argument for custom initialization of Material." % key)
                    return
        elif "material_type" in kwargs:
            self.material_type = kwargs['material_type']
            self.massdensity = config.DICT_massdensity[self.material_type]
            self.cH = config.DICT_atomic_composition[self.material_type][0]
            self.cC = config.DICT_atomic_composition[self.material_type][1]
            self.cN = config.DICT_atomic_composition[self.material_type][2]
            self.cO = config.DICT_atomic_composition[self.material_type][3]
            self.cP = config.DICT_atomic_composition[self.material_type][4]
            self.cS = config.DICT_atomic_composition[self.material_type][5]
            self.cAu = config.DICT_atomic_composition[self.material_type][6]
        else:
            logger.error("No valid arguments for Material initialization.")
            return

    def get_fX(self,element,photon_energy_eV=None):
        """
        Get the scattering factor for an element through linear interpolation.
        """
        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        SF_X = config.DICT_scattering_factors[element]
        e = constants.e
        c = constants.c
        h = constants.h
        f1 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,1])
        f2 = numpy.interp(photon_energy_eV,SF_X[:,0],SF_X[:,2])
        return complex(f1,f2) 
 
    def get_n(self,photon_energy_eV=None):
        """
        Obtains complex refractive index.
        Henke (1994): n = 1 - r_0/(2pi) lambda^2 sum_q rho_q f_q(0)
        r_0: classical electron radius
        rho_q: atomic number density of atom species q
        f_q(0): atomic scattering factor (forward scattering) of atom species q
        """

        r_0 = constants.value("classical electron radius")
        h   =  constants.h
        c   =  constants.c
        qe   = constants.e

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe

        f = self.get_f(photon_energy_eV)
        atom_density = self.get_atom_density()
        
        n = 1 - r_0/2/numpy.pi * photon_wavelength**2 * f * atom_density

        return n

    def get_dn(self,photon_energy=None):
        return (1-self.get_n(photon_energy))

    # convenience functions
    # n = 1 - delta - i beta
    def get_delta(self,photon_energy_eV=None):
        return (1-self.get_n(photon_energy_eV).real)
    def get_beta(self,photon_energy_eV=None):
        return (-self.get_n(photon_energy_eV).imag)

    def get_photoabsorption_cross_section(self,photon_energy_eV=None):

        r_0 = constants.value("classical electron radius")
        h =  constants.h
        c =  constants.c
        qe = constants.e

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe
        mu = 2*r_0*photon_wavelength*self.get_f(photon_energy_eV).imag

        return mu

    def get_transmission(self,thickness,photon_energy_eV=None):

        n = self.get_n(photon_energy_eV)
        mu = self.get_photoabsorption_cross_section(photon_energy_eV)
        rho = self.get_atom_density()

        return numpy.exp(-rho*mu*thickness)

    def get_f(self,photon_energy_eV=None):

        h  = constants.h
        c  = constants.c
        qe = constants.e

        if not photon_energy_eV:
            photon_energy_eV = self._parent._parent.source.photon.get_energy("eV")
        photon_wavelength = h*c/photon_energy_eV/qe

        atomic_composition = self.get_atomic_composition_dict()

        f_sum = 0
        for element in atomic_composition.keys():
            # sum up average atom factor
            f = self.get_fX(element,photon_energy_eV)
            f_sum += atomic_composition[element] * f
        
        return f_sum


    def get_atom_density(self):
                
        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        for element in atomic_composition.keys():
            # sum up mass
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u

        number_density = self.massdensity/M
        
        return number_density


    def get_electron_density(self):

        u = constants.value("atomic mass constant")

        atomic_composition = self.get_atomic_composition_dict()

        M = 0
        Q = 0
        for element in atomic_composition.keys():
            # sum up electrons
            M += atomic_composition[element]*config.DICT_atomic_mass[element]*u
            Q += atomic_composition[element]*config.DICT_atomic_number[element]

        electron_density = Q*self.massdensity/M
        
        return electron_density
        
        
    def get_atomic_composition_dict(self):

        atomic_composition = {}
        
        for key in self.__dict__.keys():
            if key[0] == 'c':
                exec "c_tmp = self." + key
                atomic_composition[key[1:]] = c_tmp 
 
        tmp_sum = float(sum(atomic_composition.values()))
        for element in atomic_composition.keys():
            atomic_composition[element] /= tmp_sum 
        
        return atomic_composition

