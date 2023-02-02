"""
Create pseudo-3D cube from crooked 2D seismic profile (e.g. for Petrel AVO workflow).
The user can specify the number of output inlines and the distance spacing 
between the generated inlines (which are just duplicates of the input 2D data).

NOTE: Input SEG-Y coordinates should be in projected coordinate reference system (e.g. UTM).

@author: fwrnke
@date:   28-06-2022

"""
import os 
import sys
import glob
import warnings
import argparse

import numpy as np
import segyio

#%% FUNCTIONS

def parse_trace_headers(file, filter_empty_fields=True):
    """
    Parse (non-empty) SEG-Y trace headers into dictionary.

    Parameters
    ----------
    file : segyio.segy.SegyFile
        Open file handler.
    filter_empty_fields : bool, optional
        Return only trace header fields that contain information.

    Returns
    -------
    h : dict
        Extracted trace header information.

    """
    headers = segyio.tracefield.keys
    
    h = dict()
    for k, v in headers.items():
        values = file.attributes(v)[:]
        if filter_empty_fields and not values.any():
            continue
        h[k] = values
    return h


def round_nearest(value, a):
    """ Round to nearest multiple of ``a`` """
    return round(value / a) * a


def get_azimuth(coords, kind:str = 'deg'):
    """
    Calculate dominant orientation from coordinate array (i.e. azimuth of fitted line).

    Parameters
    ----------
    coords : np.ndarray
        2D coordinate array with X and Y columns.
    kind : str, optional
        Return azimuth angle in degree (default) or radians.

    Returns
    -------
    angle : float
        Azimuth angle of fitted line through coordinates.

    """
    angle = np.arctan((coords[-1, 0] - coords[0, 0]) / (coords[-1, 1] - coords[0, 1]))
    if kind.lower() in ['deg', 'degree', 'degrees']:
        angle = np.rad2deg(angle)
    
    return angle


def rotate_coords(pts, origin=(0, 0), angle: float = 0, kind: str = 'deg'):
    """
    Rotate coordinates around origin using angle.
    Positive values rotate coordinates counter-clockwise and vis versa.

    Parameters
    ----------
    pts : np.ndarray
        2D coordinate array of points to rotate.
    origin : tuple, optional
        Rotation origin (default: (0,0)).
    angle : float, optional
        Rotation angle.
    kind : TYPE, optional
        Input format of rotation angle (default: 'deg').

    Returns
    -------
    np.ndarray
        Rotated input coordinates

    """
    if kind.lower() in ['deg', 'degree', 'degrees']:
        angle = np.deg2rad(angle)
    # rotation matrix
    R = np.array([[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle),  np.cos(angle)]])
    o = np.atleast_2d(origin)
    p = np.atleast_2d(pts)
    return np.squeeze((R @ (p.T - o.T) + o.T).T)


def wrap_text(txt, width=80):
    """
    Format textual SEG-Y header for pretty printing.
    """
    return '\n'.join([txt[i:i+width] for i in range(0, len(txt), width)])


def _isascii(txt):
    """ Check for ASCII """
    try:
        txt.decode("ascii")
    except UnicodeDecodeError:
        return False
    else:
        return True
    
def get_textual_header(path):
    """
    Read SEG-Y file and return textual header as string
    """
    with open(path, 'rb') as file:
        file.seek(0, 0) # find byte position zero relative to start of file (0)
        text_byte = file.read(3200)

    if _isascii(text_byte):
        text = text_byte.decode("ascii")
        text = wrap_text(text)
    else:
        text = text_byte.decode("cp500")  # EBCDIC encoding
        text = wrap_text(text)
    return text

def write_textual_header(path, txt:str, **kwargs_segy) -> None:
    """
    Write textual header string to SEG-Y file.
    
    Parameters
    ----------
    path : str
        SEG-Y file path.
    txt : str
        Header text.
    kwargs_segy : dict
        Optional arguments for SEG-Y I/O.

    """
    if isinstance(txt, str):
        txt = ''.join([t[:80] for t in txt.split('\n')])
    elif isinstance(txt, list):
        txt = ''.join([t[:80] for t in txt]) # silent truncating of each line if too long!
    else:
        raise ValueError(f'Not supported textual header type: {type(txt)}')
    
    header = bytes(txt, 'utf-8')
    assert len(header) == 3200, 'Binary string is too long, something went wrong...'
    
    kwargs_segy['ignore_geometry'] = True
    with segyio.open(path, 'r+', **kwargs_segy) as f:
        f.text[0] = header


def open_line(path, 
              kwargs_segyio:dict = None, 
              verbose:bool = False,
              ):
    """
    Open input 2D SEG-Y file and extract essential information.

    Parameters
    ----------
    path : str
        Input SEG-Y file path.
    kwargs_segyio : dict, optional
        Optional arguments for SEG-Y I/O (default: None).
    verbose : bool, optional
        DESCRIPTION (default: False).

    Returns
    -------
    data : np.ndarray
        Seismic data (n_traces, n_samples).
    header_traces : dict
        Dictionary of SEG-Y trace header values.
    header_binary : segyio.field.Field
        Binary header values.
    cdp_x : np.ndarray
        CDP x (easting) coordinates.
    cdp_y : TYPE
        CDP y (northing) coordinates.
    metadata : dict
        Dictionary of useful SEG-Y metadata (e.g. sampling rate).

    """
    if kwargs_segyio is None:
        kwargs_segyio = dict(strict=False, ignore_geometry=True)
    
    if verbose: print(f'[INFO]    Read SEG-Y file < {os.path.basename(path)} >')
    with segyio.open(path, mode='r' , **kwargs_segyio) as src:
        n_traces = src.tracecount           # total number of traces   
        dt = segyio.tools.dt(src) / 1000    # sample rate [ms]
        n_samples = src.samples.size        # total number of samples
        
        metadata = {'n_traces': n_traces, 'n_samples': n_samples, 'dt': dt}
        
        if verbose:
            print(f'[INFO]      n_traces:  {n_traces}')
            print(f'[INFO]      n_samples: {n_samples}')
            print(f'[INFO]      dt:        {dt} ms')
        
        # get coordinates
        cdp_x = src.attributes(segyio.TraceField.CDP_X)[:]
        cdp_y = src.attributes(segyio.TraceField.CDP_Y)[:]
        
        # get seismic data [amplitude]
        data = src.trace.raw[:]  # eager version (completely read into memory)
        
        # get trace header
        header_traces = parse_trace_headers(src)
        
        # get binary header
        header_binary = src.bin
        
    return data, header_traces, header_binary, cdp_x, cdp_y, metadata
    

def setup_cube_geometry(x, 
                        y, 
                        n_ilines:int = None,
                        spacing:float = None, 
                        origin:str = 'coordinates',
                        return_rotated_coords:bool = False, 
                        verbose:bool = False):
    """
    Create pseudo-3D cube geometry from crooked 2D line coordinates.

    Parameters
    ----------
    x : np.ndarray
        Easting (x) coordinate array.
    y : np.ndarray
        Northing (y) coordinate array.
    n_ilines : int, optional
        Number of inlines to create. Should be odd number and will be increased by 1 if not odd.
        If None: ``n_ilines`` will be calculated based on line crookedness and bin size of 2D line.
    spacing : float, optional
        Inline spacing (in CRS units, e.g. meters). If None: ``spacing`` will be calculated 
        based on line crookedness and  ``n_ilines``.
    origin : str, optional
        Rotation center coordinates. Can be 'coordinates' (default, using input coordinates) or 
        'fitted' (using fitted line coordinates).
    return_rotated_coords : bool, optional
        Return rotated cube (3D) and line (2D) coordinates for comparison (default: False).
    verbose : bool, optional
        Print additional information to STDOUT (default: False).

    Returns
    -------
    coords_cube : np.ndarray
        Array of 3D cube iline/xline coordinates with shape (n_traces x n_ilines, 2).
    coords_cube_rot : np.ndarray, optional
        Array of rotated 3D cube iline/xline coordinates with shape (n_traces x n_ilines, 2).
    coords_rot : np.ndarray
        Array of rotated input 2D coordinates with shape (n_traces, 2).
    n_ilines : int
        Number of created inlines.

    """
    x = np.asarray(x)
    y = np.asarray(y)
    assert x.size == y.size
    y_increasing = y[0] < y[-1]
    n_traces = y.size
    
    # init coordinates array
    coords = np.vstack((x, y)).T
    
    # linear fit of 2D line coordinates
    m, b = np.polyfit(x, y, 1)
    line_y = b + (m * x)
    
    # init rotation center (origin)
    if origin == 'fitted':
        origin = (x[x.size//2], line_y[x.size//2])
    elif origin == 'coordinates':
        origin = coords[coords.shape[0]//2, :]
    else:
        raise ValueError('Origin should be either ``fitted`` or ``coordinates``.')
        
    # get bin sizes from coordinate locations
    bin_sizes = np.sqrt(((coords[1:] - coords[:-1]) ** 2).sum(axis=-1))
    # select mean bin size
    bin_size = round_nearest(bin_sizes.mean(), 0.5)
    
    # get rotation of survey line (in degrees)
    azimuth = get_azimuth(coords)
    if verbose: print(f'[INFO]    Dominant profile orientation:  {azimuth:.1f} degree')
    
    # rotated 2D coordinates
    coords_rot = rotate_coords(coords, origin=origin, angle=azimuth, kind='degree')
    
    # maximum x-distance of rotated coordinates from origin
    diff_max = np.max(np.abs((coords_rot[:,0].max() - origin[0], 
                              coords_rot[:,0].min() - origin[0])))
    if n_ilines is None and spacing is None:
        n_ilines = int(np.ceil(diff_max / bin_size))
    
    if n_ilines % 2 == 1:
        n_ilines = n_ilines  
    else:
        print(f'[INFO]    Even number (> {n_ilines} <) of inlines --> increased by one.')
        n_ilines += 1
    
    if spacing is None:
        # get iline spacing from crookedness of 2D line
        spacing = np.ceil(diff_max * 2) / (n_ilines - 1)
    
    if verbose: 
        print(f'[INFO]    # ilines:        {n_ilines}')
        print(f'[INFO]    iline spacing:   {spacing:.2f} (units in CRS)')
    
    # setup cube coordinates
    xcoords = np.linspace(origin[0] - int(np.floor(n_ilines / 2) * spacing), 
                          origin[0] + int(np.floor(n_ilines / 2) * spacing), 
                          n_ilines, 
                          endpoint=True).astype(x.dtype)
    ycoords = np.linspace(coords_rot[:, 1].min() if y_increasing else coords_rot[:, 1].max(), 
                          coords_rot[:, 1].max() if y_increasing else coords_rot[:, 1].min(), 
                          n_traces, 
                          endpoint=True).astype(y.dtype)
    
    # sanity check: full coverage of crooked 2D track?
    _check_xmin = (
        np.greater(coords_rot[:,0], xcoords[0]).all() or np.less(coords_rot[:,0], xcoords[0]).all()
        )
    _check_xmax = (
        np.greater(coords_rot[:,0], xcoords[-1]).all() or np.less(coords_rot[:,0], xcoords[-1]).all()
        )
    if not all((_check_xmin, _check_xmax)):
        spacing_old = spacing
        diff_max = np.max(np.abs((coords_rot[:,0].max() - origin[0], 
                                  coords_rot[:,0].min() - origin[0])))
        spacing = np.ceil(diff_max * 2) / (n_ilines - 1)
        
        xcoords = np.linspace(origin[0] - int(np.floor(n_ilines / 2) * spacing), 
                              origin[0] + int(np.floor(n_ilines / 2) * spacing), 
                              n_ilines, 
                              endpoint=True).astype(x.dtype)
        
        print(
            f'[WARNING]    >>> 3D cube geometry based on provided spacing < {spacing_old} >',
             'does NOT cover full extend of crooked 2D line!', 
            f'It was automatically adjusted to < {spacing} > (units in CRS) <<<'
            )
    
    # setup 2D coordinate grids
    xx, yy = np.meshgrid(xcoords, ycoords)
    coords_cube_rot = np.column_stack((xx.ravel(order='F'), yy.ravel(order='F'))) 
    
    # get rotated coordinate grid for pseudo-3D cube
    coords_cube = rotate_coords(coords_cube_rot, origin=origin, angle=-azimuth)
    
    if return_rotated_coords:
        return (coords_cube, coords_cube_rot, coords_rot), n_ilines
    
    return coords_cube, n_ilines


def create_pseudo3D(fname, 
                    data, 
                    x, 
                    y, 
                    coord_scalar:int = 1, 
                    iline:int = 189, 
                    xline:int = 193, 
                    format:int = 1, 
                    dt:int = 4000, 
                    ilines=None, 
                    xlines=None, 
                    htraces:dict = None, 
                    hbinary:segyio.field.Field = None,
                    dtype:str = 'float32', 
                    overwrite:bool = True, 
                    verbose:bool = False, 
                    ) -> None:
    """
    Create and save pseudo-3D SEG-Y cube based on input 2D line.

    Parameters
    ----------
    fname : str
        Path of output SEG-Y file.
    data : np.ndarray
        Output 3D seismic data (iline, xline, samples).
    x : np.ndarray
        Cube x coordinates with shape: (n_ilines x n_traces, 2).
    y : np.ndarray
        Cube y coordinates with shape: (n_ilines x n_traces, 2).
    coord_scalar : int, optional
        Scaling factor for coordinates (default: 1).
    iline : int, optional
        Inline byte position (default: 189).
    xline : int, optional
        Crossline byte position (default: 193).
    format : int, optional
        SEG-Y format (default: 1).
    dt : int, optional
        Sampling rate in microseconds (default: 4000).
    ilines : np.ndarray, optional
        Inline numbering, must be iterable of same length as iline dimension of data. 
        If not specified, using default numbering starting with 1 (default: None).
    xlines : np.ndarray, optional
        Crossline numbering, must be iterable of same length as xline dimension of data.
        If not specified, using default numbering starting with 1 (default: None).
    htraces : dict, optional
        Dictionary of SEG-Y trace header values (default: None).
    hbinary : segyio.field.Field, optional
        Binary header values from input 2D SEG-Y file (default: None).
    dtype : str, optional
        Data output dtype (default: 'float32').
    overwrite : bool, optional
        Overwrite existing output file (default: True).
    verbose : bool, optional
        Print additional information (default: False).

    """
    # init SEG-Y headers
    HEADERS = segyio.tracefield.keys
    
    dt = int(dt)
        
    # SEG-Y data dimensions
    dims = len(data.shape)
    if dims != 3:
        raise ValueError(f'Expected 3 dimensions but {dims} were provided.')    
    nil, nxl, ns = data.shape
    
    # setup coordinates
    x = np.asanyarray(np.around(x), dtype='int32')
    y = np.asanyarray(np.around(y), dtype='int32')
    if (x.size != nil * nxl) or (y.size != nil * nxl):
        raise ValueError(
            f'Missmatch between coordinate array sizes (x: {x}, y: {y})',
            'and iline/xline number ({nil*nxl}!'
            )
    if x.size != y.size:
        raise ValueError('`x` and `y` must have same size')
    
    # init SEG-Y file specifications
    spec = segyio.spec()
    spec.iline   = iline
    spec.xline   = xline
    spec.format  = format
    spec.sorting = segyio.TraceSortingFormat.INLINE_SORTING
    
    _ilines = ilines if ilines is not None and len(ilines) == nil else list(range(1, nil + 1))
    _xlines = xlines if xlines is not None and len(xlines) == nxl else list(range(1, nxl + 1))
    spec.ilines  = _ilines
    spec.xlines  = _xlines
    spec.samples = list(range(ns))
    
    if coord_scalar is None:
        coord_scalar = 0
    if coord_scalar > 0:
        coord_scalar_mult = 1 / abs(coord_scalar)
    elif coord_scalar < 0:
        coord_scalar_mult = coord_scalar
    else:
        coord_scalar_mult = 1
    
    # prepare and filter provided trace header data from input 2D line
    if htraces is not None:
        HTRACES_NOCOPY = ['CDP_X', 'CDP_Y', 'INLINE_3D', 'CROSSLINE_3D']
        htraces = dict([(int(HEADERS.get(k)), v.astype('int32')) 
                        for k, v in htraces.items() if k not in HTRACES_NOCOPY])
    
    if os.path.isfile(fname):
        if overwrite:
            warnings.warn('File already exists and will be overwritten (overwrite=True)!')
            os.remove(fname)
        else:
            sys.exit('File already exists and overwriting is NOT permitted! Exit script.')
        
    # create SEG-Y file
    with segyio.create(fname, spec) as segyf:
        for i, il in enumerate(_ilines):
            if _ilines[0] == 0:
                il0, iln = il * nxl, (il + 1) * nxl
            else:
                il0, iln = (il - 1) * nxl, il * nxl
                            
            # update using original trace headers
            if htraces is not None:
                if verbose: print('[INFO]    Adding trace header information from input 2D line')
                segyf.header[il0:iln] = [{k: v[i] for k, v in htraces.items()} 
                                         for i in range(nxl)]
            
            # add new header values
            if verbose: print('[INFO]    Setting trace header for output SEG-Y')
            segyf.header[il0:iln] = [{
                    segyio.su.offset: 1,
                    iline: _ilines[i],
                    xline: xln,
                    segyio.su.cdpx: int(cdpx * coord_scalar_mult),
                    segyio.su.cdpy: int(cdpy * coord_scalar_mult),
                    segyio.su.ns: ns,
                    segyio.su.scalco: int(coord_scalar),
                    } 
                for xln, cdpx, cdpy in zip(_xlines, x[il0:iln], y[il0:iln])
                ]
            segyf.trace[il0:iln] = data[i, :, :].astype(dtype)
            
            if verbose and (i % 1 == 0):
                print(f'[INFO]    Writing inline < {il} >')
                # print(f'[INFO]      i: {i}, il: {il}, il0: {il0}, iln: {iln}')
        
        if hbinary is not None:
            if verbose: print('[INFO]    Adding information from input 2D line')
            segyf.bin = hbinary
            
        if verbose: print('[INFO]    Setting binary header for output SEG-Y')
        segyf.bin.update(
            tsort=segyio.TraceSortingFormat.INLINE_SORTING,
            hdt=dt, 
            dto=dt, 
            hns=ns, 
            )


def define_input_args():
    parser = argparse.ArgumentParser(
        description='Create pseudo-3D cube from crooked 2D SEG-Y file.')
    parser.add_argument('input_path', type=str, help='Input SEG-Y path.')
    parser.add_argument('--output_dir', '-o', type=str, 
                        help='Output directory for pseudo-3D SEG-Y cube(s).')
    parser.add_argument('--suffix', '-s', type=str, 
                        help='File suffix. Only used when ``input_path`` is a directory.')
    parser.add_argument('--filename_suffix', '-fns', type=str, 
                        help='Filename suffix for guided selection (e.g. "env" or "despk"). \
                            Only used when ``input_path`` is a directory.')
    parser.add_argument('--txt_suffix', type=str, default='_pseudo-3D', 
                        help='Suffix to append to output filename.')
    parser.add_argument('--n_ilines', type=int, 
                        help='Number of inlines for pseudo-3D cube (should be odd).')
    parser.add_argument('--spacing', type=int, 
                        help='Inline spacing (in units of CRS, e.g. meter).')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output SEG-Y if existing.')
    parser.add_argument('--verbose', '-V', action='store_true',
                        help='Print additional information to stdout.')
    return parser


def expand_single_segy(path_input, args, verbose) -> None:
    """ Main function """
    basepath, filename = os.path.split(path_input)
    basename, suffix = os.path.splitext(filename)
    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        output_dir = basepath
    
    if suffix.lower() not in ['.sgy', '.segy']:
        suffix = '.sgy'
    
    path_pseudo3D = os.path.join(output_dir, f'{basename}{args.txt_suffix}{suffix}')
    
    # [2D] read input SEG-Y file
    if verbose: print('[INFO] (1) Reading input 2D SEG-Y file')
    data, header_traces, header_binary, cdp_x, cdp_y, meta_segy = open_line(
        path_input, 
        verbose=verbose
        )       

    # [CUBE] setup 3D geometry
    if verbose: print('[INFO] (2) Setup output 3D geometry')
    (coords_cube, coords_cube_rot, coords_rot), n_ilines = setup_cube_geometry(
        cdp_x, cdp_y, n_ilines=args.n_ilines, spacing=args.spacing, 
        return_rotated_coords=True, verbose=verbose
        )
    
    # [CUBE] write pseudo-3D cube
    if verbose: print('[INFO] (3) Creating pseudo-3D cube')
    
    # expand data from 2D to 3D
    ## segyio: iline, xline, samples
    if verbose: print('[INFO]    Expand 2D data to 3D')
    data_3D = np.repeat(data[None, ...], n_ilines, axis=0)
    
    # write pseudo-3D cube to disk
    if verbose: print('[INFO]    Write pseudo-3D cube to disk')
    create_pseudo3D(path_pseudo3D, 
                    data_3D, 
                    x=coords_cube[:,0], 
                    y=coords_cube[:,1],
                    htraces=header_traces,
                    hbinary=header_binary, 
                    dt=meta_segy['dt'] * 1000,
                    overwrite=args.overwrite, 
                    verbose=verbose)
    
    # copy textual header from 2D to 3D SEG-Y file
    if verbose: print('[INFO]    Copy textual header')
    txt_header = get_textual_header(path_input)
    write_textual_header(path_pseudo3D, txt_header)
    

def main(input_args=None):
    parser = define_input_args()
    args = parser.parse_args(input_args)
    
    verbose = args.verbose
    
    # sanity checks
    path_input = args.input_path
    basepath, filename = os.path.split(path_input)
    basename, suffix = os.path.splitext(filename)
    
    if args.output_dir is not None:
        if not os.path.isdir(args.output_dir):
            raise IOError('Output directory does not exist.' + \
                          ' Please create directory manually before running this script!')
    
    # (1) single input file
    if os.path.isfile(path_input):
        if verbose: print(f'[INFO]    Processing file < {filename} >')
        expand_single_segy(path_input, args, verbose)
        sys.exit()
        
    # (2) input directory (multiple files)
    elif os.path.isdir(path_input):
        pattern = '*'
        pattern += f'{args.filename_suffix}' if args.filename_suffix is not None else pattern
        pattern += f'.{args.suffix}' if args.suffix is not None else '.sgy'
        file_list = sorted(glob.glob(os.path.join(path_input, pattern)))
        
        # compute files from options (2) or (3)
        if len(file_list) > 0:
            if verbose: print(f'[INFO]    Processing total of < {len(file_list)} > files')
            for file_path in file_list:
                expand_single_segy(file_path, args, verbose)
            if verbose: print(f'[SUCCESS]    Created a total of < {len(file_list)} > pseudo-3D cube(s)')
        else:
            sys.exit('No input files to process. Exit process.')
    
#%% MAIN
if __name__ == '__main__':
    main()
