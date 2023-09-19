import pytest


def pytest_addoption(parser):
    parser.addoption("--cores")


@pytest.fixture(scope='session', autouse=True)
def cores(request):
    return request.config.getoption("--cores")


class MockFastqRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality
        self.comment = ""

    def __repr__(self) -> str:
        return "|".join([self.name, self.sequence, "+", self.quality])


class MockFastaRecord:
    """Testing class used to supply a pysam FastqProxy like object"""

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self) -> str:
        return f">{self.name}\n{self.sequence}\n"
