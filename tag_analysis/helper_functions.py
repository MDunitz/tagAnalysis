import glob
import os 
import subprocess
import re
import tempfile



def _execute_r_script(r_script, success_message="R script executed successfully", verbose=False):
    if verbose:
        print("Generated R script:")
        print("="*50)
        print(r_script)
        print("="*50)

    """Helper function to execute R scripts"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        temp_r_file = f.name
    
    try:
        subprocess.run(['Rscript', temp_r_file], check=True)
        print(success_message)
    finally:
        os.unlink(temp_r_file)

def generate_file_paths_and_samples(path_to_fastq):
    """Generate all file paths and sample names in one place"""
    forward_reads = glob.glob(f"{path_to_fastq}/*R1_001*.fastq.gz")
    reverse_reads = glob.glob(f"{path_to_fastq}/*R2_001*.fastq.gz")

    forward_reads.sort()
    reverse_reads.sort()

    filtered_forward_reads = [re.sub(r"_001\.fastq", "_filtered.fastq", f) for f in forward_reads]
    filtered_reverse_reads = [re.sub(r"_001\.fastq", "_filtered.fastq", f) for f in reverse_reads]

    # Extract sample names once
    samples = []
    for f in filtered_forward_reads:
        basename = os.path.basename(f)
        sample = re.sub(r"_S[0-9]*_L.*$", "", basename)
        samples.append(sample)

    ## This is just to ensure the lists are formatted correctly for r scripts
    r_forward_reads = ', '.join([f'"{f}"' for f in forward_reads])
    r_reverse_reads = ', '.join([f'"{f}"' for f in reverse_reads])
    r_filtered_forward_reads = ', '.join([f'"{f}"' for f in filtered_forward_reads])
    r_filtered_reverse_reads = ', '.join([f'"{f}"' for f in filtered_reverse_reads])
    r_samples = ', '.join([f'"{s}"' for s in samples])
    return {
        'forward_reads': forward_reads,
        'reverse_reads': reverse_reads,
        'filtered_forward_reads': filtered_forward_reads,
        'filtered_reverse_reads': filtered_reverse_reads,
        'samples': samples,
        'r_forward_reads': r_forward_reads,
        'r_reverse_reads': r_reverse_reads,
        'r_filtered_forward_reads': r_filtered_forward_reads,
        'r_filtered_reverse_reads': r_filtered_reverse_reads,
        'r_samples': r_samples
    }



