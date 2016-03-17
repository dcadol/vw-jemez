"""
Utilities for interacting with dflow and casimir models via the virtual
watershed.
"""
import copy
import json

from numpy import fromstring, reshape
from pandas import Series


def casimir(vegetation_map, zone_map, shear_map, shear_resistance_dict):
    """
    Simple version of the CASiMiR model for vegetation succession. Before the
    model is run, we check that all the unique values from vegetation_map are
    present in the shear_resistance_dict. Otherwise the process will fail
    wherever the vegetation map value is not present in the dictionary
    on lookup.

    Arguments:
        vegetation_map (str or ESRIAsc): location on disk or ESRIAsc
            representation of the vegetation map
        shear_map (str or ESRIAsc): location on disk or ESRIAsc representation
            of the shear stress map
        shear_resistance_dict (str or dict): location on disk or dictionary
            representation of the resistance dictionary that maps
            vegetation type to shear resistance.

    Returns:
        (ESRIAsc) vegetation map updated with new values corresponding to
            succession rules
    """
    if type(vegetation_map) is str:
        vegetation_map = ESRIAsc(vegetation_map)
    elif not isinstance(vegetation_map, ESRIAsc):
        raise TypeError('vegetation_map must be type str or ESRIAsc')

    if type(shear_map) is str:
        shear_map = ESRIAsc(shear_map)
    elif not isinstance(shear_map, ESRIAsc):
        raise TypeError('shear_map must be type str or ESRIAsc')

    if type(zone_map) is str:
        zone_map = ESRIAsc(zone_map)
    elif not isinstance(zone_map, ESRIAsc):
        raise TypeError('zone_map must be type str of ESRIAsc')

    if type(shear_resistance_dict) is str:
        try:
            shear_resistance_dict = json.load(open(shear_resistance_dict))
        except ValueError:
            raise ValueError(
                'The shear_resistance_dict file is not valid JSON!'
            )
    elif not isinstance(shear_resistance_dict, dict):
        raise TypeError('shear_resistance_dict must be type str or dict')

    # init the vegetation map that will be returned
    ret_veg_map = copy.deepcopy(vegetation_map)

    for idx in range(len(shear_map.data)):
        # determine whether or not the vegetation should be reset to age zero
        shear_val_int = str(int(vegetation_map.data[idx]))

        is_not_nodata = shear_map.data[idx] != shear_map.NODATA_value

        veg_needs_reset = (
            is_not_nodata and
            shear_map.data[idx] > shear_resistance_dict[shear_val_int]
        )

        if veg_needs_reset:
            # reset vegetation to age zero while retaining veg type
            ret_veg_map.data[idx] = zone_map.data[idx]

        # whether or not the vegetation was destroyed, age by one
        if is_not_nodata:

            ret_veg_map.data[idx] += 1

    return ret_veg_map


def shear_mesh_to_asc(dflow_out_nc_path, west_easting_val, n_eastings,
                      south_northing_val, n_northings, cellsize):
    """
    Extract flow element values and locations from the dflow output netcdf
    and project these onto the grid defined by the corner of the
    grid defined by the lower-left corner of the bounding box and the
    cell size. The results are saved in ESRI .asc format to asc_out_path.

    Arguments:
        dflow_out_nc_path (str): location of the dflow netcdf on disk
        west_easting_val (float): lower-left corner easting
        south_northing_val (float): lower-left corner northing
        cellsize (float): resolution of the output grid. probably determined
            by the input vegetation map .asc

    Returns:
        (ESRIAsc) representation of gridded representation of the mesh shear
            stress data output from dflow
    """
    dflow_ds = open_dataset(dflow_out_nc_path)

    # the mesh locations are the x and y centers of the Flow (Finite) Elements
    mesh_x = dflow_ds.FlowElem_xcc
    mesh_y = dflow_ds.FlowElem_ycc
    mesh_shear = dflow_ds.taus[-1]  # take the last timestep

    x = array([west_easting_val + (i*cellsize) for i in range(n_eastings)])

    y = array([south_northing_val + (i*cellsize) for i in range(n_northings)])

    grid_x, grid_y = meshgrid(x, y)

    # use linear interp so we don't have to install natgrid
    asc_mat = griddata((mesh_x, mesh_y), mesh_shear, (grid_x, grid_y))

    # not sure why, but this makes it align with the original vegetation map
    asc_mat = flipud(asc_mat)

    data = Series(reshape(asc_mat, (n_eastings * n_northings)))

    return ESRIAsc(ncols=n_eastings, nrows=n_northings,
                   xllcorner=west_easting_val, yllcorner=south_northing_val,
                   cellsize=cellsize, data=data)


class ESRIAsc:

    def __init__(self, file_path=None, ncols=None, nrows=None,
                 xllcorner=None, yllcorner=None, cellsize=1,
                 NODATA_value=-9999, data=None):

        self.file_path = file_path
        self.ncols = ncols
        self.nrows = nrows
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.cellsize = cellsize
        self.NODATA_value = NODATA_value
        self.data = data

        # if a file is provided, the file metadata will overwrite any
        # user-provided kwargs
        if file_path:

            def getnextval(f):
                return f.readline().strip().split()[1]

            f = open(file_path, 'r')

            self.ncols = int(getnextval(f))
            self.nrows = int(getnextval(f))
            self.xllcorner = float(getnextval(f))
            self.yllcorner = float(getnextval(f))
            self.cellsize = int(getnextval(f))
            self.NODATA_value = float(getnextval(f))

            # should not be necessary for well-formed ESRI files, but
            # seems to be for CASiMiR
            data_str = ' '.join([l.strip() for l in f.readlines()])

            self.data = Series(fromstring(data_str, dtype=float, sep=' '))

            colrow_prod = self.nrows*self.ncols

            assert len(self.data) == colrow_prod, \
                "length of .asc data does not equal product of ncols * nrows" \
                "\nncols: {}, nrows: {}, ncols*nrows: {} len(data): {}".format(
                    self.ncols, self.nrows, colrow_prod, len(self.data))

    def as_matrix(self, replace_nodata_val=None):
        """
        Convenience method to give 2D numpy.ndarray representation. If
        replace_nodata_val is given, replace all NODATA_value entries with
        it.

        Arguments:
            replace_nodata_val (float): value with which to replace
                NODATA_value entries

        Returns:
            (numpy.ndarray) matrix representation of the data in the .asc
        """
        ret = copy.copy(reshape(self.data, (self.nrows, self.ncols)))
        if replace_nodata_val is not None:
            ret[ret == self.NODATA_value] = replace_nodata_val

        return ret

    def write(self, write_path):
        # replace nan with NODATA_value
        self.data = self.data.fillna(self.NODATA_value)

        with open(write_path, 'w+') as f:
            f.write("ncols {}\n".format(self.ncols))
            f.write("nrows {}\n".format(self.nrows))
            f.write("xllcorner {}\n".format(self.xllcorner))
            f.write("yllcorner {}\n".format(self.yllcorner))
            f.write("cellsize {}\n".format(self.cellsize))
            f.write("NODATA_value {}\n".format(self.NODATA_value))

            # prob not most efficient, but CASiMiR requires
            # ESRI Ascii w/ newlines
            f.write(
                '\n'.join(
                    [
                        ' '.join([str(v) for v in row])
                        for row in self.as_matrix()
                    ]
                )
            )

    def __eq__(self, other):

        if isinstance(other, ESRIAsc):
            ret = self.ncols == other.ncols
            ret = self.nrows == other.nrows and ret
            ret = self.xllcorner == other.xllcorner and ret
            ret = self.yllcorner == other.yllcorner and ret
            ret = self.cellsize == other.cellsize and ret
            ret = self.NODATA_value == other.NODATA_value and ret
            ret = all(self.data == other.data) and ret

            return ret

        return NotImplemented
