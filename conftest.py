import pytest


def pytest_addoption(parser):
    parser.addoption("--cores")


@pytest.fixture(scope='session', autouse=True)
def cores(request):
    return request.config.getoption("--cores")
