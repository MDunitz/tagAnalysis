import pandas as pd
import re
from bokeh.plotting import figure, save

from .constants import DECONTAM_PREVALENCE_THRESHOLD
from .extraction_metadata import (
    assign_extraction_batches,
    read_extraction_sheet,
    tag_id_to_batch,
)
from .helper_functions import _execute_r_script
import os

def identify_control_samples(sample_names, control_pattern=r'[Bb]lank|[Cc]ont[rol]*|[Nn]eg|DNAex|50cyc'):
    """Identify control samples based on naming pattern"""
    return [bool(re.search(control_pattern, name)) for name in sample_names]

# TODO test this function
def run_decontam_frequency(conc, counts_file_path, contam_results_file_path):
    """Run decontam using frequency method with concentrations"""
    conc_str = ', '.join(map(str, conc))
    
    r_script = f"""
    library(decontam)
    
    counts_tab <- read.csv("{counts_file_path}", row.names=1, check.names=FALSE, stringsAsFactors=FALSE,  sep='\t')
    counts_tab <- as.matrix(sapply(counts_tab, as.numeric))
    conc <- c({conc_str})

    
    # Convert to numeric matrix while preserving row and column names
    row_names <- rownames(counts_tab)
    col_names <- colnames(counts_tab)
    counts_tab <- apply(counts_tab, 2, as.numeric)
    rownames(counts_tab) <- row_names
    colnames(counts_tab) <- col_names
    
    contam_predict <- isContaminant(t(counts_tab), method='frequency', conc=conc)
    
    results <- data.frame(
        ASV = rownames(contam_predict),
        contaminant = contam_predict$contaminant,
        p_freq = contam_predict$p.freq
    )
    
    write.csv(results, "{contam_results_file_path}", row.names=FALSE)
    """
    
    _execute_r_script(r_script)
    return pd.read_csv(contam_results_file_path)

def run_decontam_prevalence(counts_file_path, contam_results_file_path, predicted_controls, threshold=0.5):
    """Run decontam using prevalence method with control samples"""
    controls_str = ', '.join(['TRUE' if x else 'FALSE' for x in predicted_controls])
    
    r_script = f"""
    library(decontam)
    
    counts_tab <- read.csv("{counts_file_path}", row.names=1, check.names=FALSE,  sep='\t')
    neg_controls <- c({controls_str})

    # Convert to numeric matrix while preserving row and column names
    row_names <- rownames(counts_tab)
    col_names <- colnames(counts_tab)
    counts_tab <- apply(counts_tab, 2, as.numeric)
    rownames(counts_tab) <- row_names
    colnames(counts_tab) <- col_names
    
    contam_predict <- isContaminant(t(counts_tab), method='prevalence', 
                                  neg=neg_controls, threshold={threshold})
    
    results <- data.frame(
        ASV = rownames(contam_predict),
        contaminant = contam_predict$contaminant,
        p_prev = contam_predict$p.prev
    )
    
    write.csv(results, "{contam_results_file_path}", row.names=FALSE)
    """
    
    _execute_r_script(r_script)
    return pd.read_csv(contam_results_file_path)

def combine_batch_p_values(batch_results, threshold):
    """Batch-aware contaminant call from per-batch prevalence results.

    decontam's batch.combine="minimum" rule: an ASV's score is its smallest
    per-batch prevalence p,

        p_ASV = min_b p_prev(ASV, b),   contaminant  <=>  p_ASV < threshold,

    so an ASV that looks like contamination in any one batch is removed.
    Batches where the ASV has no reads give p = NaN and are skipped.

    batch_results: {batch: DataFrame(ASV, p_prev)}
    Returns DataFrame(ASV, contaminant, p_prev, p_prev_batch_<b>...).
    """
    wide = pd.concat(
        {f"p_prev_batch_{b}": df.set_index("ASV")["p_prev"] for b, df in batch_results.items()},
        axis=1,
    )
    p_min = wide.min(axis=1, skipna=True)
    return pd.DataFrame(
        {"ASV": wide.index, "contaminant": (p_min < threshold).to_numpy(), "p_prev": p_min.to_numpy()}
    ).join(wide.reset_index(drop=True))


def run_decontam_prevalence_by_batch(counts_df, output_dir, predicted_controls, batches, threshold):
    """Prevalence decontam run separately within each extraction batch.

    Each batch is judged against its own kit negatives plus the shared controls
    (controls with batch None, e.g. PCR negatives, which were never extracted
    and so apply to every batch); per-batch p-values are then combined by
    `combine_batch_p_values`. A batch with no controls at all raises.
    """
    shared = [lib for lib, ctl, b in zip(counts_df.columns, predicted_controls, batches) if ctl and b is None]
    batch_results = {}
    for batch in sorted({b for b in batches if b is not None}):
        members = [lib for lib, b in zip(counts_df.columns, batches) if b == batch] + shared
        is_control = dict(zip(counts_df.columns, predicted_controls))
        controls = [is_control[lib] for lib in members]
        if not any(controls):
            raise ValueError(f"Extraction batch {batch} has no negative control and no shared control")
        subset_path = os.path.join(output_dir, f"ASVs_counts_batch_{batch}.csv")
        counts_df[members].to_csv(subset_path, sep="\t")
        batch_results[batch] = run_decontam_prevalence(
            subset_path,
            os.path.join(output_dir, f"contam_results_batch_{batch}.csv"),
            controls,
            threshold,
        )
    return combine_batch_p_values(batch_results, threshold)


def print_contaminating_taxa(contam_asvs, taxonomy_tab):
    """Print taxonomy information for contaminating ASVs"""
    if contam_asvs:
        print("Identified contaminants:")
        for asv in contam_asvs:
            if asv in taxonomy_tab.index:
                tax_string = ';'.join([str(x) for x in taxonomy_tab.loc[asv] if pd.notna(x) and str(x) != ''])
                print(f"  {asv}: {tax_string}")


# TODO consider rewriting decontam package in python to skip one round of R 
# WIll that require benchmarking/testing before publishing results that use it?
def remove_contaminants(counts_file_path, taxonomy_file_path, output_dir, 
                       conc=None, threshold=DECONTAM_PREVALENCE_THRESHOLD, clean_count_file="ASVs_counts_clean.csv",
                       extraction_sheet=None):
    """
    Remove contaminants using decontam R package
    Works much better if you have DNA conc info
    
    Parameters:
    counts_file_path: Path to the existing counts CSV file
    taxonomy_file_path: Path to the taxonomy CSV file
    output_dir: Output directory
    conc: list of concentrations (same order as samples) or None
    threshold: threshold for prevalence method
    extraction_sheet: DNA extraction sheet (CSV path or Google Sheets URL) or
        None. When given, prevalence decontam runs within each extraction
        batch against that batch's negatives (see
        run_decontam_prevalence_by_batch); libraries are placed by their TAG ID.
    """
    if extraction_sheet is not None and conc is not None:
        raise ValueError("extraction_sheet (batch-aware prevalence) and conc (frequency) are mutually exclusive")
    contam_results_file_path = os.path.join(output_dir, 'contam_results.csv')
    # Read in the data
    counts_df = pd.read_csv(counts_file_path, index_col=0, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file_path, index_col=0, sep='\t')
    relative_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    
    # Identify predicted controls
    predicted_controls = identify_control_samples(counts_df.columns)
    contam_asvs = []
    
    if any(predicted_controls) or (conc is not None and len(conc) == len(counts_df.columns)):
        
        if extraction_sheet is not None:
            print('Decontam using "prevalence" method within extraction batches')
            batches = assign_extraction_batches(
                counts_df.columns, predicted_controls,
                tag_id_to_batch(read_extraction_sheet(extraction_sheet)),
            )
            contam_results = run_decontam_prevalence_by_batch(
                counts_df, output_dir, predicted_controls, batches, threshold
            )
            contam_results.to_csv(contam_results_file_path, index=False)
        elif conc is not None and len(conc) == len(counts_df.columns):
            print('Decontam using "concentration" method')
            contam_results = run_decontam_frequency(conc, counts_file_path, contam_results_file_path)
        elif any(predicted_controls):
            print('Decontam using "prevalence" method')
            contam_results = run_decontam_prevalence(counts_file_path, contam_results_file_path, predicted_controls, threshold)
        else:
            print("No control samples found and no concentration data provided. Skipping decontamination.")
            return counts_df.copy(), relative_df.copy(), [], predicted_controls
        
        contam_asvs = contam_results[contam_results['contaminant']]['ASV'].tolist()
        
        # Print contaminating taxa
        print_contaminating_taxa(contam_asvs, taxonomy_df)
        
    else:
        print("No control samples found and no concentration data provided. Skipping decontamination.")
    
    # Remove contaminants and recalculate relative abundances
    if contam_asvs:
        print(f"Removing {len(contam_asvs)} contaminating ASVs")
        clean_counts = counts_df.drop(index=contam_asvs)
        clean_relative = clean_counts.div(clean_counts.sum(axis=0), axis=1) * 100
    else:
        print("No contaminants identified")
        clean_counts = counts_df.copy()
        clean_relative = relative_df.copy()

    # Always written (equal to the raw counts when nothing was removed): the
    # stackbar step and downstream loaders read this file unconditionally.
    clean_counts.to_csv(os.path.join(output_dir, clean_count_file), sep='\t')

    return clean_counts, clean_relative, contam_asvs, predicted_controls

