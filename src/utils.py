import os
import sys
from pathlib import Path
import glob
import subprocess
import logging
from src.config import NGSbinaries, DOCKER

logger = logging.getLogger(__name__)

def run_multiqc(output_dir, docker_config):
    '''
    '''
    cmd = ('{} run -v {}:/run_dir/ -it {}'
        ' /run_dir/. --outdir /run_dir/')\
        .format(DOCKER, docker_dict['multiqc']['image'], output_dir)
    p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    output = p1.stdout.decode('UTF-8')
    error  = p1.stderr.decode('UTF-8')

def get_fastq_files(input_dir, avoid_trimmed=False):
  '''
    Get fastq files from input directory
  '''
  fastq_list = []
  fastq_list = glob.glob(input_dir +  "/*.fastq.gz")
  if not fastq_list:
    fastq_list = glob.glob(input_dir +  "/*.fa.gz")
  if not fastq_list:
    msg = " ERROR: No input fastq files were detected"
    logging.error(msg)
    sys.exit()

  filtered_list = []
  for fastq in sorted(fastq_list):
    if avoid_trimmed == True:
        if 'trimmed' in fastq:
            continue
    if 'Undetermined' in fastq:
      continue
    else:
      filtered_list.append(fastq)

  return filtered_list

def get_bin_path(program) -> str:
    '''
        Get the PATH of a program
    '''
    path = ""
    cmd = ('which {}').format(program)
    p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    output = p1.stdout.decode('UTF-8')
    error  = p1.stderr.decode('UTF-8')
    if not error:
        if output:
            path = output.rstrip('\n')
            return path
        else:
            msg = (" ERROR: Unable to find the PATH of {}").format(program)
            logging.error(msg)
    else:
        msg = (" ERROR: Unable to find the PATH of {}").format(program)
        logging.error(msg)
