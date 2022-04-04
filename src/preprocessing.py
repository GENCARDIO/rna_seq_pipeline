import os
import sys
from pathlib import Path
import subprocess
import logging
from src.sample import Sample
from src.utils import get_fastq_files
from src.fastq import Fastq
from src.config import NGSbinaries, DOCKER
import re
logger = logging.getLogger(__name__)


def preprocess(fastq_dir, output_dir, config_dict, docker_dict):
    '''
    '''

    fastq_files = get_fastq_files(fastq_dir, avoid_trimmed=True)

    # This list will hold sample objects
    sample_list = []

    # Check fastq files are well formated
    seen_fq = set()

    for fq in sorted(fastq_files):

        f_pair = Fastq(fq, expect_paired=True, force_naming_convention=False)
        fq1 = f_pair.fq1
        fq2 = f_pair.fq2
        sample_name = f_pair.sample_name

        sample_folder = os.path.join(output_dir, sample_name)

        if sample_name in seen_fq:
            continue

        seen_fq.add(sample_name)

        # Create a new Sample object
        sample = Sample(sample_name)

        if not os.path.isdir(sample_folder):
            os.mkdir(sample_folder)

        sample.add("sample_folder", sample_folder)
        fastq_folder = os.path.join(sample_folder, "FASTQ_FOLDER")

        if not os.path.isdir(fastq_folder):
            os.mkdir(fastq_folder)

        sample.add("fastq_folder", fastq_folder)
        sample.add("sample_name", sample_name)
        trimmed_fq1, trimmed_fq2 = fastp(sample_name, fastq_folder,
            fq1, fq2, config_dict['threads'], NGSbinaries.FASTP)
        sample.add("ready_fq1", trimmed_fq1)
        sample.add("ready_fq2", trimmed_fq2)

        fastqc_report_fq1 = fastqc(fq1, fastq_folder, config_dict['threads'], docker_dict)
        fastqc_report_fq2 = fastqc(fq2, fastq_folder, config_dict['threads'], docker_dict)

        sample_list.append(sample)
    return sample_list

def fastqc(fq, output_dir, threads, docker_dict) -> str:
    '''
    '''
    cmd = ('{} run -v {}:/fastq_dir/ -v {}:/output_dir/'
        ' --rm {} -t -f fastq -o /fastq_dir/. /output_dir/{}')\
        .format(DOCKER, os.path.dirname(fq), output_dir, docker_dict['fastqc']['image'],
            threads, os.path.basename(fq))

    fastqc_report_name = fq.replace(".fastq.gz", "") + "_fastqc.zip"
    fastqc_report = os.path.join(output_dir, fastqc_report_name)

    if not os.path.isfile(fastqc_report):

        msg = (" INFO: Trimming sample {}").format(sample_name)
        logging.info(msg)
        p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')

    return fastqc_report

def fastp(sample_name, output_dir, fq1, fq2, threads, fastp_exe):
    '''
        Trim raw FASTQ files using fastp

        :param str sample_name: sample name extracted from fastq
        :param str output_dir: output directory
        :param str fq1: raw fastq1 (R1)
        :param str fq2: raw fastq2 (R2)
        :param str threads: number of cpu cores
        :param str fastp_exe: fastp binary location

        :returns: trimmed_fq1, trimmed_fq2
        :rtype: tuple
    '''

    fq1_name = os.path.basename(fq1)
    fq2_name = os.path.basename(fq2)
    trimmed_fq1_name = fq1_name.replace(".fastq.gz", ".trimmed.fastq.gz")
    trimmed_fq2_name = fq2_name.replace(".fastq.gz", ".trimmed.fastq.gz")

    trimmed_fq1 = os.path.join(output_dir, trimmed_fq1_name)
    trimmed_fq2 = os.path.join(output_dir, trimmed_fq2_name)

    #json_name = ("{}{}").format(sample_name, ".json")
    json_name = "fastp.json"
    output_json = os.path.join(output_dir, json_name)

    # Now trimming with fastp
    bashCommand = ('{} -i  {} -I {} -o {} -O {} -w {} -j {}') \
      .format(fastp_exe, fq1, fq2, trimmed_fq1, trimmed_fq2, threads, output_json)
    if not os.path.isfile(trimmed_fq1) and not os.path.isfile(trimmed_fq2):

      msg = (" INFO: Trimming sample {}").format(sample_name)
      logging.info(msg)
      logging.info(bashCommand)

      p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
      output = p1.stdout.decode('UTF-8')
      error  = p1.stderr.decode('UTF-8')
      if error:
        if re.search("error", error):
          msg = (" ERROR: Something wrong happened with fastp trimming for sample {}")\
            .format(sample_name)
          logging.error(error)
          sys.exit()
        else:
          msg = (" INFO: FASTQ Trimming ended successfully for sample {}")\
            .format(sample_name)
          logging.info(msg)
      else:
        msg = (" INFO: FASTQ Trimming ended successfully for sample {}")\
            .format(sample_name)
        logging.info(msg)
    else:
      msg = (" INFO: Skipping trimming for sample {}").format(sample_name)
      logging.info(msg)

    return trimmed_fq1, trimmed_fq2
