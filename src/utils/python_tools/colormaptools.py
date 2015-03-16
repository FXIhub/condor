#====================#
# Python tools - CXI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com


from matplotlib import colors
from pylab import *

cmaps = {}

cdict_jet = {'red':   [(0.0,0.0,0.0),
                       (0.375,0.0,0.0),
                       (0.625,1.0,1.0),
                       (0.875,1.0,1.0),
                       (1.0,0.5,0.5)],
             
             'green':  [(0.0,0.0,0.0),
                        (0.125,0.0,0.0),
                        (0.375,1.0,1.0),
                        (0.625,1.0,1.0),
                        (0.875,0.0,0.0),
                        (1.0,0.0,0.0)],
             
             'blue':  [(0.0,0.5,0.5),
                       (0.125,1.0,1.0),
                       (0.375,1.0,1.0),
                       (0.625,0.0,0.0),
                       (1.0,0.0,0.0)]}

cmaps["myjet"] = matplotlib.colors.LinearSegmentedColormap('myjet', cdict_jet, 1024)

cdict_jet_lightbg = {'red':   [(0.0,1.0,1.0),
                               (0.125,0.0,0.0),
                               (0.375,0.0,0.0),
                               (0.625,1.0,1.0),
                               (0.875,1.0,1.0),
                               (1.0,0.5,0.5)],
                     
                     'green':  [(0.0,1.0,1.0),
                                (0.125,0.0,0.0),
                                (0.375,1.0,1.0),
                                (0.625,1.0,1.0),
                                (0.875,0.0,0.0),
                                (1.0,0.0,0.0)],
                     
                     'blue':  [(0.0,1.0,1.0),
                               (0.125,1.0,1.0),
                               (0.375,1.0,1.0),
                               (0.625,0.0,0.0),
                               (1.0,0.0,0.0)]}

cmaps["jet_lightbg"] = matplotlib.colors.LinearSegmentedColormap('jet_lightbg', cdict_jet_lightbg, 1024)


cdict_jet_lightbg2 = {'red':   [(0.0,1.0,1.0),
                                (0.15,0.8,0.8),
                                (0.225,0.0,0.0),
                                (0.375,0.0,0.0),
                                (0.625,1.0,1.0),
                                (0.875,1.0,1.0),
                                (1.0,0.5,0.5)],
                     
                      'green':  [(0.0,1.0,1.0),
                                 (0.15,0.8,0.8),
                                 (0.225,0.0,0.0),
                                 (0.375,1.0,1.0),
                                 (0.625,1.0,1.0),
                                 (0.875,0.0,0.0),
                                 (1.0,0.0,0.0)],
                     
                      'blue':  [(0.0,1.0,1.0),
                                (0.15,1.0,1.0),
                                (0.225,1.0,1.0),
                                (0.375,1.0,1.0),
                                (0.625,0.0,0.0),
                                (1.0,0.0,0.0)]}

cmaps["jet_lightbg2"] = matplotlib.colors.LinearSegmentedColormap('jet_lightbg2', cdict_jet_lightbg2, 1024)

cdict_redblue = {'red':   [(0.00,0.0,0.0),
                           (0.25,1.0,1.0),
                           (0.50,1.0,1.0),
                           (0.75,0.0,0.0),
                           (1.00,0.0,0.0)],
                     
                 'green':  [(0.00,0.0,0.0),
                            (0.25,0.0,0.0),
                            (0.50,1.0,1.0),
                            (0.75,0.0,0.0),
                            (1.00,0.0,0.0)],
                 
                 'blue':  [(0.00,0.0,0.0),
                           (0.25,0.0,0.0),
                           (0.50,1.0,1.0),
                           (0.75,1.0,1.0),
                           (1.00,0.0,0.0)]}

cmaps["redblue"] = matplotlib.colors.LinearSegmentedColormap('redblue', cdict_redblue, 1024)

cdict_grayalpha = {'red':   [(0.0,0.0,0.0),
                             (1.0,0.0,0.0)],
                   
                   'green':  [(0.0,0.0,0.0),
                              (1.0,0.0,0.0)],
                   
                   'blue':  [(0.0,0.0,0.0),
                             (1.0,0.0,0.0)],
                   
                   'alpha': [(0.0,0.0,0.0),
                             (1.0,1.0,1.0)]}

cmaps["grayalpha"] = matplotlib.colors.LinearSegmentedColormap('grayalpha', cdict_grayalpha, 1024)




def make_colorbar_blank(filename,Nx,Ny,colormap=cm.jet,orientation="vertical"):
    X,Y = meshgrid(arange(Nx),arange(Ny))
    if orientation == "vertical":
        C = -Y
    else:
        C = X
    imsave(filename,C,cmap=colormap)

def make_colorbar(filename,**kwargs):
    from matplotlib import pyplot
    import matplotlib as mpl
    cmap = kwargs.get("cmap",cm.jet)
    label= kwargs.get("label","")
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)
    Nx = kwargs.get("Nx",100)
    Ny = kwargs.get("Ny",400)
    dpi= kwargs.get("dpi",400)
    scaling = kwargs.get("scaling","linear")
    orientation = kwargs.get("orientation","vertical")
    fig = pyplot.figure(figsize=(Nx/100.,Ny/100.))
    if orientation == "vertical":
        ax = fig.add_axes([0.15, 0.1, 0.15, 0.8])
    else:
        ax = fig.add_axes([0.1, 0.15, 0.8, 0.15])
    if scaling == "linear":
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if scaling == "log":
        if vmin <= 0 or vmin == None:
            vmin = .1
        if vmax <= 0 or vmax == None:
            vmax = 10
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,norm=norm,orientation=orientation)
    cb.set_label(label)
    pyplot.savefig(filename,dpi=dpi)


def complex_array_to_rgb(X, theme='dark', rmax=None):
    '''Takes an array of complex number and converts it to an array of [r, g, b],
    where phase gives hue and saturaton/value are given by the absolute value.
    Especially for use with imshow for complex plots.'''
    absmax = rmax or abs(X).max()
    Y = zeros(X.shape + (3,), dtype='float')
    Y[..., 0] = angle(X) / (2 * pi) % 1
    if theme == 'light':
        Y[..., 1] = np.clip(np.abs(X) / absmax, 0, 1)
        Y[..., 2] = 1
    elif theme == 'lightlog':
        Y[..., 1] = np.log10(np.clip(10*(np.abs(X) / absmax*0.9+0.1), 1, 10))
        Y[..., 2] = 1
    elif theme == 'dark':
        Y[..., 1] = 1
        Y[..., 2] = np.clip(np.abs(X) / absmax, 0, 1)
    elif theme == 'darklog':
        Y[..., 1] = 1
        Y[..., 2] = np.log10(np.clip(10*(np.abs(X) / absmax*0.9+0.1), 1, 10))
    Y = matplotlib.colors.hsv_to_rgb(Y)
    return Y
