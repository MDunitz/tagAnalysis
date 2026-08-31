import glob
import os 
import subprocess
from .constants import RANKS
from .helper_functions import _execute_r_script
import pandas as pd

# Extract/cleanup
def remove_primers_cutadapt(data_file_path, fwd, rev, rev_rc, fwd_rc, 
                           output_directory="fastq_cutadapt", min_length=100):
    """Remove primers using cutadapt"""
    os.makedirs(output_directory, exist_ok=True)
    data_files = glob.glob(f"{data_file_path}/*R1*")
    
    for R1 in data_files:
        R2 = R1.replace('_R1_', '_R2_')
        
        R1_basename = os.path.basename(R1)
        R2_basename = os.path.basename(R2)
        
        cmd = [
            'cutadapt',
            '-a', f'^{fwd}...{rev_rc}',
            '-A', f'^{rev}...{fwd_rc}',
            '--discard-untrimmed',
            '-m', str(min_length),
            '-o', f'{output_directory}/{R1_basename}',
            '-p', f'{output_directory}/{R2_basename}',
            R1,
            R2
        ]
    
        print(f"Processing {R1_basename}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            pass
            # print(f"Successfully processed {R1_basename}")
        else:
            print(f"Error processing {R1_basename}: {result.stderr}")

def filter_and_trim_reads(
        r_forward_reads, 
        r_reverse_reads, 
        r_filtered_forward_reads, 
        r_filtered_reverse_reads, 
        output_dir,
        maxEE=(2,2), 
        minLen=150, 
        truncLen=(230,200), 
        rm_phix=True, 
        filtered_results_file="filter_results.csv"):
    filtered_results_path = os.path.join(output_dir, filtered_results_file)
    """Filter and trim reads, return summary DataFrame"""
    
    r_script = f"""
    library(dada2)

    forward_reads <- c({r_forward_reads})
    reverse_reads <- c({r_reverse_reads})
    filtered_forward_reads <- c({r_filtered_forward_reads})
    filtered_reverse_reads <- c({r_filtered_reverse_reads})

    filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                    reverse_reads, filtered_reverse_reads, maxEE=c{maxEE},
                    rm.phix={str(rm_phix).upper()}, minLen={minLen}, truncLen=c{truncLen})

    result <- cbind(filtered_out, perc_kept = filtered_out[,2]/filtered_out[,1])
    write.csv(result, "{filtered_results_path}", row.names=TRUE)
    print(result)
    """
    
    _execute_r_script(r_script)
    return pd.read_csv(filtered_results_path, index_col=0)

def learn_errors_and_denoise(r_filtered_forward_reads, r_filtered_reverse_reads, r_samples, output_dir,
                           multithread=8, pool="pseudo", 
                           err_fwd_reads_file="err_forward_reads.rds",
                           err_rev_reads_file="err_reverse_reads.rds",
                           dada_fwd_file="dada_forward.rds",
                           dada_rev_file="dada_reverse.rds",
                           binned_quality_bins=None,
                           verbose=False,
                           ):
    """Learn error rates and perform denoising.

    binned_quality_bins: for instruments that emit binned quality scores
        (NextSeq/NovaSeq), the list of bin values Illumina collapses Q-scores
        onto, e.g. [2, 12, 24, 40]. When provided, learnErrors uses a binned
        error-estimation function (makeBinnedQualErrfun) so the error model is
        learned correctly for binned data. When None (default), the standard
        loess error function is used -- correct for unbinned MiSeq data, so
        older datasets are unaffected. Confirm the exact bins against the
        quality-profile heatmap per run; Illumina's binning can change.
    """
    verbose = str(verbose).upper()

    # Build the errorEstimationFunction argument. For binned data, learnErrors
    # must be told the bins via makeBinnedQualErrfun; otherwise the default
    # (loess) is used by omitting the argument entirely.
    if binned_quality_bins is not None:
        bins_r = ", ".join(str(int(b)) for b in binned_quality_bins)
        err_fn_setup = f"binQ <- makeBinnedQualErrfun(c({bins_r}))"
        err_fn_arg = ", errorEstimationFunction=binQ"
    else:
        err_fn_setup = ""
        err_fn_arg = ""

    r_script = f"""
    library(dada2)

    filtered_forward_reads <- c({r_filtered_forward_reads})
    filtered_reverse_reads <- c({r_filtered_reverse_reads})
    samples <- c({r_samples})

    {err_fn_setup}

    # Learn error rates
    err_forward_reads <- learnErrors(filtered_forward_reads, multithread={multithread}{err_fn_arg})
    err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread={multithread}{err_fn_arg})

    # Dereplicate
    derep_forward <- derepFastq(filtered_forward_reads, verbose={verbose})
    derep_reverse <- derepFastq(filtered_reverse_reads, verbose={verbose})

    names(derep_forward) <- samples
    names(derep_reverse) <- samples

    # DADA2 denoising step
    dada_forward <- dada(derep_forward, err=err_forward_reads, pool="{pool}", multithread={multithread})
    dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="{pool}", multithread={multithread})

    # Save results
    saveRDS(err_forward_reads, "{os.path.join(output_dir, err_fwd_reads_file)}")
    saveRDS(err_reverse_reads, "{os.path.join(output_dir, err_rev_reads_file)}")
    saveRDS(dada_forward, "{os.path.join(output_dir, dada_fwd_file)}")
    saveRDS(dada_reverse, "{os.path.join(output_dir, dada_rev_file)}")
    """
    
    _execute_r_script(r_script)


# Transform data
def merge_pairs_and_remove_chimeras(
        r_filtered_forward_reads, 
        r_filtered_reverse_reads, 
        r_samples,
        output_file_dir,
        trimOverhang=True, 
        minOverlap=12,  
        dada_fwd_file="dada_forward.rds",
        dada_rev_file="dada_reverse.rds", 
        seq_file = "sequence_table.csv",
        seq_no_chimeras_file ="sequence_table_nochim.csv",
        merged_fraction_file="merged_fractions.csv",
        non_chimeric_fraction_file="nonchimeric_fraction.csv",
        rds_seq_file="seqtab_nochim.rds",
        rds_amplicon_file="merged_amplicons.rds",
        verbose=False
        ):
    dada_fwd_file_path = os.path.join(output_file_dir, dada_fwd_file)
    dada_rev_file_path = os.path.join(output_file_dir, dada_rev_file)
    """Merge paired reads and remove chimeras"""
    verbose = str(verbose).upper()
    r_script = f"""
    library(dada2)

    # Load previously saved objects
    dada_forward <- readRDS("{dada_fwd_file_path}")
    dada_reverse <- readRDS("{dada_rev_file_path}")

    filtered_forward_reads <- c({r_filtered_forward_reads})
    filtered_reverse_reads <- c({r_filtered_reverse_reads})
    samples <- c({r_samples})

    # Recreate dereplicated objects (needed for mergePairs)
    derep_forward <- derepFastq(filtered_forward_reads, verbose={verbose})
    derep_reverse <- derepFastq(filtered_reverse_reads, verbose={verbose})
    names(derep_forward) <- samples
    names(derep_reverse) <- samples

    # Merge paired reads
    merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, 
                                   trimOverhang={str(trimOverhang).upper()}, minOverlap={minOverlap})

    # Make sequence table
    seqtab <- makeSequenceTable(merged_amplicons)

    print("Merged fraction:")
    merged_fractions <- rowSums(seqtab) / rowSums(makeSequenceTable(dada_forward))
    print(merged_fractions)

    # Remove chimeras
    # TODO pull out method?
    seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = T, method = 'consensus')

    print("Non-chimeric fraction:")
    nonchim_fraction <- sum(seqtab.nochim)/sum(seqtab)
    print(nonchim_fraction)

    # Save results
    write.csv(seqtab, "{os.path.join(output_file_dir, seq_file)}")
    write.csv(seqtab.nochim, "{os.path.join(output_file_dir, seq_no_chimeras_file)}")
    write.csv(merged_fractions, "{os.path.join(output_file_dir, merged_fraction_file)}")
    write.csv(nonchim_fraction, "{os.path.join(output_file_dir, non_chimeric_fraction_file)}")
    
    saveRDS(merged_amplicons, "{os.path.join(output_file_dir, rds_amplicon_file)}")
    saveRDS(seqtab.nochim, "{os.path.join(output_file_dir, rds_seq_file)}")
    """
    
    _execute_r_script(r_script)
    
    # Return key metrics
    merged_fractions = pd.read_csv(os.path.join(output_file_dir, merged_fraction_file), index_col=0)
    nonchim_fraction = pd.read_csv(os.path.join(output_file_dir, non_chimeric_fraction_file), index_col=0)
    
    return {
        'merged_fractions': merged_fractions,
        'nonchim_fraction': nonchim_fraction.iloc[0,0]
    }

def create_summary_table(
        r_samples, 
        output_file_dir,
        dada_fwd_file="dada_forward.rds", 
        dada_rev_file="dada_reverse.rds", 
        rds_seq_table="seqtab_nochim.rds", 
        filtered_file="filter_results.csv", 
        rds_amplicon_file="merged_amplicons.rds", 
        read_count_file="read-count-tracking.tsv",
        verbose=False):
    verbose = str(verbose).upper()
    dada_fwd_file_path = os.path.join(output_file_dir, dada_fwd_file)
    dada_rev_file_path = os.path.join(output_file_dir, dada_rev_file)
    rds_seq_table_file_path = os.path.join(output_file_dir, rds_seq_table)
    filtered_file_file_path = os.path.join(output_file_dir, filtered_file)
    rds_amplicon_file_file_path = os.path.join(output_file_dir, rds_amplicon_file)
    read_count_file_path = os.path.join(output_file_dir, read_count_file)

    """Create summary table using R script - only uses R for DADA2-specific operations"""
    
    r_script = f"""
    library(dada2)

    # Load saved objects
    dada_forward <- readRDS("{dada_fwd_file_path}")
    dada_reverse <- readRDS("{dada_rev_file_path}")
    seqtab.nochim <- readRDS("{rds_seq_table_file_path}")
    merged_amplicons <- readRDS("{rds_amplicon_file_file_path}")
    filtered_out <- read.csv("{filtered_file_file_path}", row.names=1)
    
    samples <- c({r_samples})

    # Define getN function and create summary table
    getN <- function(x) sum(getUniques(x))

    summary_tab <- data.frame(row.names=samples, samp=samples, input=filtered_out[,1],
                   filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                   dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                   nonchim=rowSums(seqtab.nochim),
                   final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

    write.table(summary_tab, "{read_count_file_path}", quote=FALSE, sep="\t", col.names=NA)
    print(summary_tab)
    """
    
    _execute_r_script(r_script)
    return pd.read_csv(read_count_file_path, sep="\t", index_col=0)

## TODO pull out code below into separate file

def prepare_plotting_data(output_file_dir, summary_tab_file="read-count-tracking.tsv"):
    # Read the summary table
    summary_tab = pd.read_csv(os.path.join(output_file_dir, summary_tab_file), sep="\t", index_col=0)
    
    # Melt the data (exclude the final percentage column)
    columns_to_melt = [col for col in summary_tab.columns if col != 'final_perc_reads_retained']
    melted_data = pd.melt(
        summary_tab[columns_to_melt].reset_index(),
        id_vars=['samp'],
        var_name='Step',
        value_name='Reads'
    )
    melted_data['Reads'] = pd.to_numeric(melted_data['Reads'], errors='coerce')
    melted_data = melted_data.dropna()

    # Prepare percentage data
    pct_data = melted_data.copy()
    
    # Calculate percentages relative to input for each sample
    for sample in pct_data['samp'].unique():
        sample_mask = pct_data['samp'] == sample
        input_reads = pct_data[(sample_mask) & (pct_data['Step'] == 'input')]['Reads'].iloc[0]
        pct_data.loc[sample_mask, 'Reads'] = pct_data.loc[sample_mask, 'Reads'] / input_reads
    
    # Add sample type classification
    pct_data['guess_type'] = 'sample'
    control_pattern = r'[Bb]lank|[Cc]ont[rol]*|[Nn]eg|DNAex|50cyc'
    control_mask = pct_data['samp'].str.contains(control_pattern, regex=True, na=False)
    pct_data.loc[control_mask, 'guess_type'] = 'control'
    
    return melted_data, pct_data

def create_asv_outputs(output_dir, seqtab_nochim_file="sequence_table_nochim.csv", asv_map_file="asv_mapping.csv", asvs_fasta_file="ASVs.fa"):
    seqtab_nochim_file_path = os.path.join(output_dir, seqtab_nochim_file)
    counts_file_path = os.path.join(output_dir, "ASVs_counts.csv")
    asvs_fa_file_path = os.path.join(output_dir, asvs_fasta_file)

    seqtab_nochim = pd.read_csv(seqtab_nochim_file_path, index_col=0)
    
    asv_seqs = seqtab_nochim.columns.tolist()
    asv_headers = [f"ASV_{i+1}" for i in range(len(asv_seqs))]
    
     # Create mapping file
    mapping_df = pd.DataFrame({
        'ASV_ID': asv_headers,
        'sequence': asv_seqs
    })
    mapping_path = os.path.join(output_dir, asv_map_file)
    mapping_df.to_csv(mapping_path, sep='\t', index=False)

    with open(asvs_fa_file_path, 'w') as f:
        for header, seq in zip(asv_headers, asv_seqs):
            f.write(f">{header}\n{seq}\n")
    
    # Fast transpose and save
    asv_tab = seqtab_nochim.T
    asv_tab.index = asv_headers
    asv_tab.to_csv(counts_file_path, sep='\t')
    
    return mapping_df, asv_tab

# threshold=40 means 40% confidence as some lineages would be unclassified (default 60%)
def assign_taxonomy(output_dir, mapping_df, reference_db_path, taxonomy_file="ASV_taxonomy.csv", ranks=RANKS, threshold=40, processors=16):

    taxonomy_results_file_path = os.path.join(output_dir, taxonomy_file)

    asv_ids = mapping_df['ASV_ID'].tolist()
    sequences = mapping_df['sequence'].tolist()

    r_sequences = ', '.join([f'"{seq}"' for seq in sequences])
    r_asv_ids = ', '.join([f'"{asv_id}"' for asv_id in asv_ids])

    """Create summary table using R script - only uses R for DADA2-specific operations"""
    
    r_script = f"""
    library(DECIPHER)
    library(dada2)
    load("{reference_db_path}")
    
    sequences <- c({r_sequences})
    asv_ids <- c({r_asv_ids})

    # Stupid hack to use the asv_id mapped sequences instead of using r getSequences function
    sequences <- gsub("U", "T", sequences)  # Convert U to T
    sequences <- gsub("[^ATGC]", "N", sequences)  # Replace invalid chars with N
    
    
    # Create DNAStringSet with clean names
    dna <- DNAStringSet(sequences)
    names(dna) <- asv_ids

    
    # Classify
    tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="top", threshold={threshold}, processors={processors})

    # Create results with clean ASV IDs
    results <- data.frame(
        ASV_ID = asv_ids,
        sequence = sequences,
        taxonomy = sapply(tax_info, function(x) paste(x$taxon[-1], collapse=";")),
        confidence = sapply(tax_info, function(x) min(x$confidence)),
        stringsAsFactors = FALSE
    )
    
    write.table(results, "{taxonomy_results_file_path}", sep='\t', row.names=FALSE, quote=FALSE)
    print(paste("Classified", length(asv_ids), "ASVs"))

    """
    
    _execute_r_script(r_script)
    taxonomy_df = pd.read_csv(taxonomy_results_file_path, sep="\t")
    # Apply function to create rank columns
    taxonomy_df[ranks] = taxonomy_df.apply(split_taxonomy_to_ranks, axis=1)
    
    # Set ASV_ID as index for final output
    taxonomy_df.set_index('ASV_ID', inplace=True)
    taxonomy_df.to_csv(taxonomy_results_file_path, sep="\t", na_rep='')
    return taxonomy_df

# todo is there a cleaner way to rewrite this? 
def split_taxonomy_to_ranks(row, ranks=RANKS):
        """Split taxonomy string into rank columns"""
        taxonomy_string = row['taxonomy']
        
        if pd.notna(taxonomy_string) and taxonomy_string.strip():
            # Split taxonomy string
            taxa_parts = taxonomy_string.split(';')
            
            # Pad with last taxon if needed
            if len(taxa_parts) > 0:
                last_taxon = taxa_parts[-1]
                while len(taxa_parts) < len(ranks):
                    taxa_parts.append(last_taxon)
            
            # Take only first len(ranks) elements
            taxa_parts = taxa_parts[:len(ranks)]
            
            # Create result series with rank names as index
            result = pd.Series(index=ranks, dtype=str)
            for i, rank in enumerate(ranks):
                if i < len(taxa_parts) and taxa_parts[i].strip():
                    result[rank] = taxa_parts[i].strip()
                else:
                    result[rank] = ""
            
            return result
        else:
            # Fill with empty strings if no taxonomy
            return pd.Series([""] * len(ranks), index=ranks)
 
def prepare_data_for_contamination_plot(relative_df, contam_asvs, predicted_controls):
        # Prepare data for plotting
    plot_data = []
    
    for sample in relative_df.columns:
        library_type = "Control" if predicted_controls[list(relative_df.columns).index(sample)] else "Real"
        
        for asv in relative_df.index:
            
            plot_data.append({
                'Sample': sample,
                'ASV': asv,
                'relabund': relative_df.loc[asv, sample],
                'decontam': asv in contam_asvs,
                'library_type': library_type
            })

    plot_df = pd.DataFrame(plot_data)
    # Aggregate by contamination status and sample
    agg_data = plot_df.groupby(['Sample', 'decontam', 'library_type'])['relabund'].sum().reset_index()
    
    return agg_data


def prepare_relative_abundance_data(counts_file_path, taxonomy_file_path):
    """
    Transform counts and taxonomy data into long format for plotting
    
    Parameters:
    counts_file_path: Path to counts CSV file
    taxonomy_file_path: Path to taxonomy CSV file
    
    Returns:
    relative_long_df: DataFrame in long format with relative abundances and taxonomy
    """
    
    # Read data
    counts_df = pd.read_csv(counts_file_path, index_col=0, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file_path, index_col=0, sep='\t')
    
    # Calculate relative abundances
    relative_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    
    # Melt to long format (equivalent to R's melt function)
    relative_long_df = relative_df.reset_index().melt(
        id_vars='index', 
        var_name='sample', 
        value_name='relabund'
    ).rename(columns={'index': 'ASV'})
    
    # Add taxonomy information (equivalent to R's cbind with tax_df)
    taxonomy_cols = [col for col in taxonomy_df.columns if col not in ['sequence', 'taxonomy', 'confidence']]
    relative_long_df = relative_long_df.merge(
        taxonomy_df[taxonomy_cols].reset_index(), 
        left_on='ASV', 
        right_on='ASV_ID', 
        how='left'
    )
    
    return relative_long_df