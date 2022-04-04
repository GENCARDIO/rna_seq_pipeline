import os
import sys
import argparse
import logging
import re
from pathlib import Path
import logging
from src.preprocessing import preprocess
from src.config import load_genome_config, load_docker_config
from src.map import align_reads
from src.quantification import apply_featureCounts

main_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(main_dir+"/src")

def parse_args():
    '''
    '''
    parser = argparse.ArgumentParser(description="Run a basic rna-seq pipeline")
    parser.add_argument("--fastq_dir",
        type=str, required=True, help="Input fastq file directory")
    parser.add_argument("--output_dir", help="Output directory",
        type=str, required=True)
    parser.add_argument("--config_yaml", help="Config yaml",
        type=str, required=True)
    parser.add_argument("--docker_yaml", help="Docker images yaml",
        type=str, required=True)
    parser.add_argument("-t", "--threads", type=int, default=4,
        help="Num. of CPU threads to operate", dest='threads')
    parser.add_argument("-r", "--reference", required=True, type=str,
        choices=['hg19', 'hg38'])

    return parser.parse_args()

if __name__== "__main__":

    logging.getLogger().setLevel(logging.INFO)
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    args = parse_args()
    fastq_dir  = args.fastq_dir
    output_dir = args.output_dir

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    config_dict = load_genome_config(args.config_yaml)
    docker_dict = load_docker_config(args.docker_yaml)

    args_dict = vars(args)
    config_dict = {**config_dict, **args_dict}

    # Fastq preprocessing
    sample_list = preprocess(fastq_dir, output_dir, config_dict, docker_dict)

    sample_list = align_reads(sample_list, config_dict, docker_dict)

    sample_list = apply_featureCounts(sample_list, config_dict, docker_dict)
