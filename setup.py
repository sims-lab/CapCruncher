from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(

rust_extensions=[RustExtension("capcruncher.libcapcruncher","Cargo.toml", debug=False, binding=Binding.PyO3)],

)
