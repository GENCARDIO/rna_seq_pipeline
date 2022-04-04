# rna_seq_pipeline

## Graphical description

![Graphical description](/img/rnaseq_pipeline.png)

## Installation
This code has been successfully tested to run in python 3.8.
You can check your python3 version by prompting: `python3 -V`.
Also, you will need `pip3` as a package manager.

It is advisable to use a virtual environment to build the project:
`python3 -m virtualenv mypython`

### 1. Docker

Docker is required tto install and run some of the pipeline components.
To intall it you can follow these instructions:
(https://docs.docker.com/engine/install/)

### 2. Install python3 dependencies
`pip3 install -r requeriments.txt`

## Basic usage:

First, change the configuration parameters from `config.yaml` that handle genome files location

The, you can run the pipeline with the following command
`python3 rna_seq_pipeline.py --threads <NUM_CPUS> --fastq_dir <FASTQ_DIR> --output_dir <OUTPUT_DIR> -r <hg19/hg38> --config_yaml config.yaml --docker_yaml docker_images.yaml`
