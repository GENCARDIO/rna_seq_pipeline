import os
import sys
from pathlib import Path
import subprocess
import logging
from src.sample import Sample
from src.utils import get_fastq_files
from src.fastq import Fastq
import re
from src.config import NGSbinaries

logger = logging.getLogger(__name__)

class InvalidBAM(Exception):
    pass


def salmon_alignment(sample_list, config_dict, docker_dict):
    '''
    '''
    for sample in sample_list:
        pass



def index_bam(bam_in) -> None:
    '''
    '''
    bai = bam_in + ".bai"
    cmd = "{} index {}".format(NGSbinaries.SAMTOOLS, bam_in)

    if not os.path.isfile(bai):
        p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        if error:
            raise InvalidBAM(error)

def mark_duplicates(bam_in) -> str:
    '''
        Markduplicates with Picard
        :param:
    '''

    bam_out = bam_in.replace(".bam", ".rmdup.bam")
    picard_metrics = bam_out.replace(".bam", ".picard.txt")

    cmd = "java -jar {} MarkDuplicates -I {} -O {} -M {}".format(
        NGSbinaries.PICARD, bam_in, bam_out,picard_metrics)

    if not os.path.isfile(bam_out):
        msg = (" INFO: Marking duplicates for {}").format(bam_in)
        logging.info(msg)

        p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
    else:
        msg = (" INFO: Skipping duplicate marking for {}").format(bam_in)
        logging.info(msg)

    index_bam(bam_out)

    return bam_out

def align_reads(sample_list, config_dict, docker_dict):
    '''
    '''
    for sample in sample_list:

        bam_folder = os.path.join(sample.sample_folder, "BAM_FOLDER")
        if not os.path.isdir(bam_folder):
            os.mkdir(bam_folder)
        sample.add("bam_folder", bam_folder)

        hisat2 = Hisat2(sample.name, sample.ready_fq1, sample.ready_fq2,
            config_dict['GRCh38']['hisat2_index'], bam_folder)
        bam = hisat2.align()
        sample.add("raw_bam", bam)

        rmdup_bam = mark_duplicates(bam)
        sample.add("ready_bam", rmdup_bam)

    return sample_list

class Hisat2():
    '''
    '''
    def __init__(self, sample_name, fq1, fq2, genome_index, output_dir, threads=2):
        self._sample_name = sample_name
        self._fq1 = fq1
        self._fq2 = fq2
        self._genome_index = genome_index
        self._output_dir = output_dir
        self._threads = threads
        self._bam = output_dir +"/"+ self._sample_name + ".bam"

    def align(self) -> str:
        '''
            Setting --rna-strandness RF for strand-specific library
        '''

        read_group = "--rg-id={} --rg SM:{} --rg PL:ILLUMINA".format(self._sample_name,
            self._sample_name)
        summary_file = self._bam.replace(".bam", ".summary.alignment.txt")

        cmd = ('{} -x {} -1 {} -2 {} -p {} {} --summary-file {} --rna-strandness RF'
            ' | {} view -Sb - | {} sort -T {} -o {}')\
            .format(NGSbinaries.HISAT2, self._genome_index, self._fq1 ,
            self._fq2, self._threads, read_group, summary_file, NGSbinaries.SAMTOOLS,
            NGSbinaries.SAMTOOLS, self._sample_name, self._bam)

        if not os.path.isfile(self._bam):
            msg = (" INFO: Trimming sample {}").format(self._sample_name)
            logging.info(msg)

            p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')
        else:
            msg = (" INFO: Skipping hisat2 mapping for sample {}").format(self._sample_name)
            logging.info(msg)

        self.index_bam(self._bam)
        return self._bam

    @staticmethod
    def index_bam(bam_in) -> None:
        '''
        '''

        bai = bam_in + ".bai"
        cmd = "{} index {}".format(NGSbinaries.SAMTOOLS, bam_in)

        if not os.path.isfile(bai):
            p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')
            if error:
                raise InvalidBAM(error)

class STAR():
    '''
    '''
    def __init__(self, sample_name, fq1, fq2, genome_index, output_dir):
        self._sample_name = sample_name
        self._fq1 = fq1
        self._fq2 = fq2
        self._genome_index = genome_index
        self._output_dir = output_dir
        self._bam = output_dir +"/"+ self._sample_name + ".bam"
