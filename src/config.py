import os
import sys
import yaml
from yaml import SafeLoader
import logging
from collections import defaultdict
import subprocess
from src.utils import get_bin_path
logger = logging.getLogger(__name__)

main_dir = os.path.dirname(os.path.abspath(__file__))
binaries_dir = os.path.join(main_dir, "../binaries")
DOCKER = get_bin_path("docker")

class NGSbinaries():
    HISAT2   = os.path.join(binaries_dir, "hisat2-2.2.1-Linux_x86_64/hisat2-2.2.1/hisat2")
    SAMTOOLS = os.path.join(binaries_dir, "samtools/bin/samtools")
    FASTP    = os.path.join(binaries_dir, "fastp")
    PICARD   = os.path.join(binaries_dir, "picard.jar")

def load_docker_config(docker_yaml) -> dict():
    '''
    '''
    with open(docker_yaml) as f:
        docker_dict = yaml.load(f, Loader=SafeLoader)

    for program in docker_dict:
        image = docker_dict[program]['image']
        validate_image(program, image)
    return docker_dict

def validate_image( program, image):
    '''
    '''

    DOCKER = get_bin_path("docker")
    cmd = '{} image ls {}'.format(DOCKER, image)
    p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    output = p1.stdout.decode('UTF-8')
    error  = p1.stderr.decode('UTF-8')
    if not error:
        if output:
            c_lines = 0
            for line in output.split('\n'):
                c_lines +=1
            if c_lines !=3:
                msg = (" ERROR: docker image {} was not found").format(image)
                logging.error(msg)
            else:
                msg = (" INFO: found docker image {}").format(image)
                logging.info(msg)
        else:
            msg = (" ERROR: docker image {} was not found").format(image)
            logging.error(msg)


def load_genome_config(config_yaml = None) -> dict():
    '''
    '''
    with open(config_yaml) as f:
        config_dict = yaml.load(f, Loader=SafeLoader)
    return config_dict

# class Config():
#     def __init__(self, config_yaml = None):
#
#         # Check that specified config file exist
#         if config_yaml is None:
#             msg = " ERROR: missing configuration yaml"
#             logging.error(msg)
#             raise FileNotFoundError(msg)
#
#         Config._CONFIG_YAML= config_yaml
#
#         with open(Config._CONFIG_YAML) as f:
#             Config._CONFIG = yaml.load(f, Loader=SafeLoader)

    # def get_property(self, property_name):
    #     '''
    #     '''
    #     if property_name not in Config._CONFIG:
    #         return None
    #
    #     return Config._CONFIG[property_name]

        # convert dict keys to class attributes
        # if Config._CONFIG is not None:
        #     for a, b in Config._CONFIG.items():
        #         if isinstance(b, (list, tuple)):
        #            setattr(self, str(a), [Config(x) if isinstance(x, dict) else x for x in b])
        #         else:
        #            setattr(self, str(a), Config(b) if isinstance(b, dict) else b)
