import pandas as pd
import re
from bokeh.plotting import figure, save

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
                       conc=None, threshold=0.5, clean_count_file="ASVs_counts_clean.csv"):
    """
    Remove contaminants using decontam R package
    Works much better if you have DNA conc info
    
    Parameters:
    counts_file_path: Path to the existing counts CSV file
    taxonomy_file_path: Path to the taxonomy CSV file
    output_dir: Output directory
    conc: list of concentrations (same order as samples) or None
    threshold: threshold for prevalence method
    """
    contam_results_file_path = os.path.join(output_dir, 'contam_results.csv')
    # Read in the data
    counts_df = pd.read_csv(counts_file_path, index_col=0, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_file_path, index_col=0, sep='\t')
    relative_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    
    # Identify predicted controls
    predicted_controls = identify_control_samples(counts_df.columns)
    contam_asvs = []
    
    if any(predicted_controls) or (conc is not None and len(conc) == len(counts_df.columns)):
        
        if conc is not None and len(conc) == len(counts_df.columns):
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
        
        # Save cleaned data
        clean_counts_file = os.path.join(output_dir, clean_count_file)
        clean_counts.to_csv(clean_counts_file, sep='\t')
    else:
        print("No contaminants identified")
        clean_counts = counts_df.copy()
        clean_relative = relative_df.copy()
    
    return clean_counts, clean_relative, contam_asvs, predicted_controls

