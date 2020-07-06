from setuptools import setup, find_packages

 setup(
   name='ccanalyser',
   version='0.0.1',
   author='asmith, dsims',
   author_email='alastair.smith@ndcls.ox.ac.uk',
   packages=find_packages(),
   entry_points={'console_scripts': ['ccanalyser = ccanalyser.capturec.cli:main',
                                     'cc_pipeline = ccanalyser.capturec.capturec_pipeline:run_pipeline']
                },
   include_package_data=True,
   url='https://github.com/sims-lab/capture-c.git',
   license='LICENSE',
   description='Performs complete processing of capture-c data',
   long_description=open('README.txt').read(),
   python_requires='>=3.6',
)
