# import os
# from pybedtools.bedtool import BedTool
# import pytest
# import glob
# import click


# # Pre-run setup
# dir_test = os.path.realpath(os.path.dirname(__file__))
# dir_test_run = os.path.join(dir_test, 'pipeline_test_run')
# dir_package = os.path.dirname(dir_test)
# dir_data = os.path.join(dir_package, "data")

# def test_is_on():

#     from capcruncher.utils import is_on
#     assert is_on('True') == True
#     assert is_on('F') == False

# def test_is_none():

#     from capcruncher.utils import is_none
#     assert is_none('True') == False
#     assert is_none('X') == False
#     assert is_none('') == True
#     assert is_none('none') == True
#     assert is_none(None) == True

# def test_get_human_readable_number_of_bp():

#     from capcruncher.utils import get_human_readable_number_of_bp
#     assert get_human_readable_number_of_bp(999) == '999bp'
#     assert get_human_readable_number_of_bp(1000) == '1.0kb'

# def test_get_re_site():

#     from capcruncher.utils import get_re_site

#     assert get_re_site('GATC') == 'GATC'
#     assert get_re_site('dpnii') == 'GATC'

#     with pytest.raises(ValueError):
#         get_re_site('XXXXXX')

# def test_format_coordinates():
#     from capcruncher.utils import format_coordinates

#     bed = format_coordinates('chr1:1000-2000')
#     assert isinstance(bed, BedTool)
#     assert len(bed) == 1

#     with pytest.raises(ValueError):
#         bed = format_coordinates('chrXXXXX1:1000-2000')
    
#     with pytest.raises(ValueError):
#         bed = format_coordinates('chr1:1000:2000')

#     test_slices = os.path.join(dir_data, "test", "test_slices.bed")
#     bed = format_coordinates(test_slices)
#     assert isinstance(bed, BedTool)
#     assert len(bed) == 4

#     with pytest.raises(ValueError):
#         test_bed_bad = os.path.join(dir_data, "test", "test_capture_bad_format.bed")
#         format_coordinates(test_bed_bad)
