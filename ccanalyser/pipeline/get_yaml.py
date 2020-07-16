import os
import sys
import shutil


SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)

shutil.copy(f'{PACKAGE_DIR}/data/capturec_pipeline.yml', 'capturec_pipeline.yml')
