import requests
import json
import subprocess
import os


import subprocess
import os

def run_boltz_prediction(input_path, use_msa=True, out_dir="./boltz_results"):
    """
    Executes the Boltz predict CLI command from Python.
    """

    command = ["boltz", "predict", input_path]
    
    if use_msa:
        command.append("--use_msa_server")
    
    command.extend(["--out_dir", out_dir])

    # 2. Execute the command
    try:
        print(f"Running: {' '.join(command)}")
        result = subprocess.run(
            command, 
            capture_output=True, 
            text=True, 
            check=True
        )
        
        return {
            "status": "success",
            "stdout": result.stdout,
            "output_directory": os.path.abspath(out_dir),
            "command_executed": " ".join(command)
        }

    except subprocess.CalledProcessError as e:
        return {
            "status": "error",
            "message": e.stderr,
            "return_code": e.returncode
        }


import yaml

def update_yaml_sequence(input_file, new_sequence, output_file):
    with open(input_file, 'r') as f:
        config = yaml.safe_load(f)

    config['input'][0]['entity'][0]['protein'] = new_sequence

    # 3. Save it back
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return output_file