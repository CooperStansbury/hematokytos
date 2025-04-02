import os
import sys
import argparse
import logging
import traceback
from pathlib import Path
import yaml
import json
import shutil
import subprocess

def run_command(command, capture_output=False, text=True):
    """
    Runs a command and streams the output to the console in real-time.

    Args:
        command (str or list): The command to run as a string or a list of arguments.
        capture_output (bool): Whether to capture the command's output. Set to False for real-time output.
        text (bool): Whether to decode the output as text.

    Returns:
        int: The return code of the command.
    """
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=text,
            shell=True,
        )

        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())

        return process.returncode

    except Exception as e:
        print(f"Command failed with error: {e}")
        return 1  # Return a non-zero exit code to indicate failure


def create_output_directory(config, overwrite):
    """
    Creates the output directory specified in the config.

    Args:
        config (dict): A dictionary containing configuration settings, including 'paths' and 'output_directory'.
        overwrite (bool): If True, overwrite an existing directory. If False, raise an error if the directory exists.
    """
    output_dir = config['paths']['output_directory']

    if os.path.exists(output_dir):
        if overwrite:
            try:
                shutil.rmtree(output_dir)
                logging.info(f"Overwriting existing output directory: {output_dir}")
            except Exception as e:
                logging.error(f"Failed to overwrite existing output directory: {e}")
                raise
        else:
            logging.error(f"Output directory already exists: {output_dir}")
            raise ValueError(f"Output directory '{output_dir}' already exists. Use --overwrite to replace.")

    try:
        working_dir = f"{output_dir}work/"
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(working_dir, exist_ok=True)
        logging.info(f"Created output directory: {output_dir}")
        logging.info(f"Created working directory: {working_dir}")
        return working_dir
    except Exception as e:
        logging.error(f"Failed to create output directory: {e}")
        raise


def load_yaml(filepath):
    """
    Loads a YAML file and returns its content as a Python dictionary.
    """
    try:
        with open(filepath, 'r') as file:
            data = yaml.safe_load(file)
        return data
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        return None


def setup_logging():
    """Sets up logging to both console and file."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),  # Log to console
        ]
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the Epi2me wf-single-cell workflow.")
    parser.add_argument(
        "--config_path",
        dest="config_path",
        default="./config.yaml",
        help="Path to the input configurtion file",
    )

    parser.add_argument(
        "--overwrite",
        dest="force",
        action="store_true", # Make --force a boolean flag
        help="If present, clear and overwrite existing output directory",
    )

    args = parser.parse_args()
    config_path = args.config_path
    force = args.force
    setup_logging()

    logging.info(f"Config file path: {config_path}")
    logging.info(f"Force overwrite: {force}")

    """ Load the config """
    config = load_yaml(config_path)
    logging.info("Loaded configuration:")
    logging.info(json.dumps(config, indent=4)) 

    """ Set up the directories """
    working_dir = create_output_directory(config, force)

    """ Set up the moldules """
    return_code = run_command("module load openjdk")
    logging.info(f"Loaded openjdk")
    return_code = run_command("module load singularity")
    logging.info(f"Loaded singularity")
    nextflow_path = config['paths']['nextflow']
    
    return_code = run_command(f"export PATH='{nextflow_path}:$PATH'")
    logging.info(f"Set Nextflow path: {nextflow_path}")

    """ set up the nextflow command """
    fastq_fpaths = config['paths']['fastq_paths']
    fastq_paths = [x.strip() for x in open(fastq_fpaths) if not x.startswith("#")]
    logging.info(f"Running pipeline for {len(fastq_paths)} FASTQs from: {fastq_fpaths}")

    fastq_paths_str = " ".join(fastq_paths)

    # data paths
    output_dir = config['paths']['output_directory']
    ref_dir = config['paths']['ref_genome_dir']

    # nextflow params
    process_name = config['nextflow']['name']
    flow_config = config['nextflow']['nextflow_config']
    numba_config = config['nextflow']['numba_config']

    # pipeline params
    params = []
    for k, v in config['pipeline'].items():
        params.append(f" --{k} {v}")    
    param_string = " ".join(params)

    nextflow_command = (
        f"nextflow run epi2me-labs/wf-single-cell"
        f" --fastq {fastq_paths_str}"  
        f" --ref_genome_dir {ref_dir}"
        f" --out_dir {output_dir}"
        f" -w {working_dir}"
        f" -c {numba_config}"
        f" -name {process_name}"
        f" -with-report"
        f" -profile singularity"
    ) + param_string

    logging.info(f" ---------- Prepared Nextflow command: ")
    logging.info(f"{nextflow_command}")
    return_code = run_command(nextflow_command)
    

    
    


