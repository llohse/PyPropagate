
def get_metric_prefix(numbers):
    import numpy as np
                  
    if not isinstance(numbers,list):
        numbers = list(numbers)
        
    def get_exponent(number):
        return int(np.log10(np.abs(number))) if number != 0 else 0
    
    from units import metric_prefixes
    
    exponents = [get_exponent(number) for number in numbers]
    largest = max(exponents,key=lambda x:abs(x))
    
    closest = min(metric_prefixes, key=lambda x:abs(x[1]+1-largest))
    return (closest[2],10**closest[1])

def get_unitless_bounds(array):

    from .units import get_unit

    bounds = []
        
    for l,r in array.bounds:
        unit = get_unit(l)
        if unit == None:
            unit = get_unit(r)
        try:
            if unit == None:
                bounds.append((float(l),float(r),1)) 
            else:
                bounds.append((float(l/unit),float(r/unit),unit)) 
        except:
            raise ValueError('Cannot convert to unitless expression: %s with unit: %s' % ((l,r),unit))
    
    return bounds
    
def image_plot(carr,ax = None,figsize = None,title = None, **kwargs):
    import matplotlib.pyplot as plt

    # fix missing \text support
    from expresso.pycas import latex as rlatex
    latex = lambda x:rlatex(x).replace(r'\text',r'\mathrm')

    fig = None
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    if title:
        ax.set_title(title)
    
    e = get_unitless_bounds(carr)
    xprefix,xfactor = get_metric_prefix(e[1][:2])
    yprefix,yfactor = get_metric_prefix(e[0][:2])
        
    extent = [float(e[1][0])/xfactor,float(e[1][1])/xfactor,float(e[0][0])/yfactor,float(e[0][1])/yfactor]

    if 'aspect' not in kwargs:
        kwargs['aspect'] = 'auto'
    if 'origin' not in kwargs:
        kwargs['origin'] = 'lower'

    image = ax.imshow(carr.data, extent=extent, **kwargs )
    ax.set_ylabel("$%s$ [$%s %s$]" % (latex(carr.axis[0]),yprefix,latex(e[0][2])))
    ax.set_xlabel("$%s$ [$%s %s$]" % (latex(carr.axis[1]),xprefix,latex(e[1][2])))
    
    if fig:
        fig.colorbar(image)

    if ax == None:
        plt.show()

    return image
        
def line_plot(carr,ax = None,ylabel = None,figsize = None,title = None,**kwargs):
    import matplotlib.pyplot as plt
    import numpy as np

    # fix missing \text support
    from expresso.pycas import latex as rlatex
    latex = lambda x:rlatex(x).replace(r'\text',r'\mathrm')

    fig = None
    if ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    
    if title:
        ax.set_title(title)
    
    e = get_unitless_bounds(carr)[0]
    
    prefix,factor = get_metric_prefix(e[:2])
    
    lines = ax.plot(np.linspace(float(e[0])/factor,float(e[1])/factor,carr.data.shape[0]),carr.data, **kwargs)
    ax.set_xlabel("$%s$ [$%s %s$]" % (latex(carr.axis[0]),prefix,latex(e[2])))
    if ylabel: ax.set_ylabel(ylabel)

    if ax == None:
        plt.show()

    return lines[0]

def expression_to_field(expression,settings):
    import expresso.pycas
    import numpy as np

    s = settings.simulation_box
    expr = settings.get_optimized(expression)
    sym = expresso.pycas.get_symbols_in(expr)

    from .coordinate_ndarray import CoordinateNDArray

    if sym - {s.x,s.y,s.z} != set():
        raise ValueError('cannot create field: contains non coordinate symbols %s' % ','.join([str(a) for a in sym - {s.x,s.y,s.z}]))
    if len(sym) == 0:
        c = complex(expr)
        if c.imag == 0:
            return c.real
        return c
        #raise ValueError('cannot create field: expression contains no symbols')
    elif len(sym) == 1:
        x = sym.pop()
        keys = tuple([getattr(s,p % x.name) for p in ['%smin','%smax','N%sd']])
        xmin,xmax,nx = settings.get_numeric( keys )
        nxmin,nxmax = settings.get_as( (xmin,xmax) , float )
        nx = settings.get_as( nx , int )
        npx = np.linspace(nxmin,nxmax,nx)
        data =  expresso.pycas.numpyfy(expr)(**{x.name:npx})
        res =  CoordinateNDArray(data,[(xmin,xmax)],(x,),settings.get_numeric_transform())
    elif len(sym) == 2:
        y,x = sorted([sym.pop(),sym.pop()],key = lambda x:x.name)[::-1]
        keys = tuple([getattr(s,p % i) for i in (x,y) for p in ['%smin','%smax','N%sd']])
        xmin,xmax,nx,ymin,ymax,ny = settings.get_numeric( keys )
        nxmin,nxmax,nymin,nymax = settings.get_as( (xmin,xmax,ymin,ymax) , float )
        nx,ny = settings.get_as( (nx,ny) , int )
        npy,npx = np.meshgrid(np.linspace(nymin,nymax,ny),np.linspace(nxmin,nxmax,nx))
        data =  expresso.pycas.numpyfy(expr)(**{x.name:npx,y.name:npy})
        res =  CoordinateNDArray(data,[(xmin,xmax),(ymin,ymax)],(x,y),settings.get_numeric_transform())
    else:
        raise ValueError('cannot create field: three dimensional field creation not implemented')

    return res

def plot(arg, *args, **kwargs):
    """
    Simple plot function for 1D and 2D coordinate arrays. If the data is complex, the absolute square value of the data will be plottted.
    
    Parameters
    -----------
    arg: coordinate array
          the input data
    
    **kwargs: additional parameters to be passed to the plot functions
    
    Returns
    --------
    plot: output of ax.plot for 1D and ax.imshow for 2D arrays
    
    """
    import expresso.pycas
    import numpy as np
    from coordinate_ndarray import CoordinateNDArray

    if isinstance(arg,expresso.pycas.Expression):
        from .settings import Settings

        if len(args) > 0 and isinstance(args[0],Settings):
            settings = args[0]
            args = list(args)
            del args[0]
        else:
            settings = kwargs.get('settings')
            if not settings:
                raise ValueError('cannot plot expression: no settings provided')
            del kwargs['settings']

        arg = expression_to_field(arg,settings)

        if isinstance(arg,(float,complex)):
            print arg
            return

    elif not isinstance(arg,CoordinateNDArray):
        raise ValueError('cannot plot non CoordinateNDArray object. For plotting regular arrays please use the matplotlib.pyplot module.')

    if not np.can_cast(arg.data.dtype, np.float128):
        if np.all(arg.data.imag == 0): arg = arg.real
        else: arg = abs(arg) ** 2
    if len(arg.axis) == 1: return line_plot(arg, *args, **kwargs)
    elif len(arg.axis) == 2: return image_plot(arg, *args, **kwargs)
    else: raise ValueError("input array must be one or two dimensional")
    
    
  