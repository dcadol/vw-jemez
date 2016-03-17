"""
Tests for posting/ingetsting dflow and casimir data to each respective model
to/from the virtual watershed.
"""
import json
import numpy
import unittest

from jemez.dflow_casimir import ESRIAsc, casimir


class TestDflow(unittest.TestCase):
    """
    Functions for working with DFLOW inputs and outputs
    """
    def setUp(self):
        self.ascii_veg = 'test/data/vegcode.asc'
        self.excel_veg_to_nval = \
            'test/data/lookup_table.xlsx'
        self.expected_ascii_roughness = \
            'test/data/roughness.asc'

        self.expected_ascii_nvals = \
            ESRIAsc(self.expected_ascii_roughness)

    def test_vegmap_properly_read(self):

        vegmap_mat = ESRIAsc(self.ascii_veg).as_matrix()

        vmat_unique = numpy.unique(vegmap_mat)
        vmat_expected = numpy.array([-9999, 100, 101, 102, 106, 210, 215],
                                    dtype='f8')

        assert (vmat_unique == vmat_expected).all()

    def test_casimir(self):

        # load the expected ESRIAsc output from running casimir
        expected_output = ESRIAsc(
            'test/data/expected_veg_output.asc'
        )

        # test results when loaded from file
        veg_map_file = self.ascii_veg
        shear_map_file = 'test/data/shear.asc'
        casimir_required_data = 'test/data/resist_manning_lookup.xlsx'
        zone_map_file = 'test/data/zonemap.asc'

        generated_output = casimir(veg_map_file, zone_map_file, shear_map_file,
                                   casimir_required_data)

        assert expected_output == generated_output, \
            "expected: {}\ngenerated: {}".format(
                expected_output.as_matrix(), generated_output.as_matrix()
            )

        # test results when using ESRIAsc instances
        veg_map = ESRIAsc(veg_map_file)
        zone_map = ESRIAsc(zone_map_file)
        shear_map = ESRIAsc(shear_map_file)

        generated_output = casimir(veg_map, zone_map,
                                   shear_map, casimir_required_data)

        assert expected_output == generated_output, \
            "expected: {}\ngenerated: {}".format(
                expected_output.as_matrix(), generated_output.as_matrix()
            )

    def test_casimir_with_dflow_io(self):
        assert False

    def test_mesh_to_asc(self):
        assert False
