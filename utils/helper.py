import os
import errno
import tarfile
import time
import yaml
import shutil
import glob
from subprocess import Popen, PIPE
import multiprocessing
from multiprocessing.pool import ThreadPool


def get_values_from_input(yamldata, key):

    values = {str(sample): yamldata[sample][key] for sample in yamldata}
    return values

