import os
from .helper_functions import generate_file_paths_and_samples
from .plotting import create_read_count_plots, generate_quality_profile_plots
from .etl import filter_and_trim_reads, learn_errors_and_denoise, merge_pairs_and_remove_chimeras, create_summary_table, prepare_plotting_data


def run_dada2_pipeline(path_to_fastq_files, path_to_output_dir, dataset_name, img_dir=None, **kwargs):
    """Main pipeline function that orchestrates all steps"""
    if img_dir is None:
        img_dir = os.path.join(path_to_output_dir, "imgs")
    
    # Step 1: Generate all file paths and sample names
    print("Step 1: Generating file paths and sample names...")
    file_info = generate_file_paths_and_samples(path_to_fastq_files)
    print(f"Found {len(file_info['samples'])} samples: {file_info['samples']}")
    
    # Step 2: Generate quality profile plots
    print("Step 2: Generating quality profile plots...")
    generate_quality_profile_plots(
        file_info['forward_reads'], 
        file_info['reverse_reads'], 
        dataset_name,
        path_to_output_dir,
        plot_count=kwargs.get('quality_plot_count', 3),
        img_dir=img_dir,
    )
    
    # TODO (OPTIONAL) quantitatively figure out trim lengths


    # Step 3: Filter and trim reads
    print("Step 3: Filtering and trimming reads...")
    filter_summary = filter_and_trim_reads(
        file_info['r_forward_reads'],
        file_info['r_reverse_reads'],
        file_info['r_filtered_forward_reads'],
        file_info['r_filtered_reverse_reads'],
        path_to_output_dir,
        maxEE=kwargs.get('maxEE', (2,2)),
        minLen=kwargs.get('minLen', 150),
        truncLen=kwargs.get('truncLen', (230,200)),
        rm_phix=kwargs.get('rm_phix', True),
    )
    print("Filter summary:")
    print(filter_summary)
    
    # Step 4: Learn errors and denoise
    print("Step 4: Learning error rates and denoising...")
    learn_errors_and_denoise(
        file_info['r_filtered_forward_reads'],
        file_info['r_filtered_reverse_reads'],
        file_info['r_samples'],
        output_dir=path_to_output_dir,
        multithread=kwargs.get('multithread', 8),
        pool=kwargs.get('pool', "pseudo"),
        binned_quality_bins=kwargs.get('binned_quality_bins', None)
    )
    
    # Step 5: Merge pairs and remove chimeras
    print("Step 5: Merging pairs and removing chimeras...")
    merge_results = merge_pairs_and_remove_chimeras(
        file_info['r_filtered_forward_reads'],
        file_info['r_filtered_reverse_reads'],
        file_info['r_samples'],
        output_file_dir=path_to_output_dir,
        trimOverhang=kwargs.get('trimOverhang', True),
        minOverlap=kwargs.get('minOverlap', 12)
    )
    print(f"Non-chimeric fraction: {merge_results['nonchim_fraction']:.3f}")
    
    # Step 6: Create summary table
    print("Step 6: Creating summary table...")
    summary_table = create_summary_table(file_info['r_samples'], output_file_dir=path_to_output_dir)
    print("\nSummary table:")
    print(summary_table)
    
    # Step 7: Prepare plotting data and create plots
    print("Step 7: Creating read retention plots...")
    melted_data, pct_data = prepare_plotting_data(output_file_dir=path_to_output_dir)
    read_count_plot_file_path = os.path.join(img_dir, kwargs.get('plot_output', "read_count_tracking.html"))
    create_read_count_plots(melted_data, pct_data, output_file=read_count_plot_file_path)
    print("\nDADA2 pipeline completed successfully!")
    
    return {
        'file_info': file_info,
        'filter_summary': filter_summary,
        'merge_results': merge_results,
        'summary_table': summary_table,
        'melted_data': melted_data,
        'pct_data': pct_data
    }



