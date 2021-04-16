import os
import pytest
import glob


# Pre-run setup
dir_test = os.path.realpath(os.path.dirname(__file__))
dir_test_run = os.path.join(dir_test, 'pipeline_test_run')
dir_package = os.path.dirname(dir_test)

def test_is_on():

    from ccanalyser.utils import is_on
    assert is_on('True') == True
    assert is_on('F') == False

def test_is_none():

    from ccanalyser.utils import is_none
    assert is_none('True') == False
    assert is_none('X') == False
    assert is_none('') == True
    assert is_none('none') == True
    assert is_none(None) == True

def test_get_human_readable_number_of_bp():

    from ccanalyser.utils import get_human_readable_number_of_bp
    assert get_human_readable_number_of_bp(999) == '999bp'
    assert get_human_readable_number_of_bp(1000) == '1.0kb'

def test_get_re_site():

    from ccanalyser.utils import get_re_site

    assert get_re_site('GATC') == 'GATC'
    assert get_re_site('dpnii') == 'GATC'

    with pytest.raises(ValueError):
        get_re_site('XXXXXX')
