import os
from functools import wraps


def write_script(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        outdir = kwargs.get("outdir")
        if outdir is None:
            raise ValueError("Error: The 'outdir' parameter is missing or None.")

        job_name = kwargs.get("job_name", "job_script")
        threads = kwargs.get("n_jobs", 1)

        script_path = os.path.join(outdir, "Script")
        os.makedirs(script_path, exist_ok=True)

        log_path = os.path.join(outdir, "Log")
        os.makedirs(log_path, exist_ok=True)

        stderr_path = os.path.join(log_path, f"{job_name}.e")
        stdout_path = os.path.join(log_path, f"{job_name}.o")

        script_header = f"""#!/bin/bash
#$ -N {job_name}
#$ -S /bin/bash
#$ -V
#$ -l hostname=Client1
#$ -e {stderr_path}
#$ -o {stdout_path}
#$ -pe smp {threads}

echo "================================================="
echo "Started on       : $(date)"
echo "Running on node  : $(hostname)"
echo "Current job ID   : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "================================================="
"""

        script_footer = """
echo "================================================="
echo "Finished on      : $(date)"
echo "================================================="
echo "Done."
"""

        job_cmd = func(*args, **kwargs)

        script_file = os.path.join(script_path, f"{job_name}.sh")
        with open(script_file, "w") as fh_script:
            fh_script.write(script_header)
            fh_script.write(job_cmd + "\n")
            fh_script.write(script_footer)

        print(f"Script written at: {script_file}")
        
        return job_cmd

    return wrapper