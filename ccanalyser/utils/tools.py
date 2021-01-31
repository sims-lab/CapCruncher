import multiprocessing as mp
import re
from collections import Counter

import pandas as pd
import numpy as np
from pysam import FastxFile
from xopen import xopen
from ccanalyser.utils.helpers import get_re_site
from random import randint
from multiprocessing import Queue, Pool, Manager, Process
from typing import Union



