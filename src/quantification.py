import os
import sys
from pathlib import Path
import subprocess
import logging
import re
from src.sample import Sample
from src.config import NGSbinaries, DOCKER
logger = logging.getLogger(__name__)


def salmon_quantification(sample_list, config_dict, docker_dict):
    '''
    '''
    pass

def apply_featureCounts(sample_list, config_dict, docker_dict):
    '''
    '''
    for sample in sample_list:

        count_file_name = sample.name + ".counts.txt"
        count_file = os.path.join(sample.bam_folder, count_file_name)
        sample.add("count_file", count_file)
        cmd = ('{} run -v {}:/gtf_dir/ -v {}:/bam_dir/ -v {}:/out_dir/'
            ' {} featureCounts -a /gtf_dir/{} -t exon '
            '-g gene_id -o /out_dir/{} /bam_dir/{}').format(
                DOCKER,
                os.path.dirname(config_dict['GRCh38']['gtf']),
                sample.bam_folder,
                sample.bam_folder,
                docker_dict['featureCounts']['image'],
                os.path.basename(config_dict['GRCh38']['gtf']),
                os.path.basename(count_file_name),
                os.path.basename(sample.ready_bam))
        if not os.path.isfile(count_file):
            p1 = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

    return sample_list
