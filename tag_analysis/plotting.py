import os
import pandas as pd
from bokeh.plotting import save, figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.transform import jitter
from bokeh.palettes import Category20
from bokeh.layouts import column
import seaborn as sns
import iqplot

from .helper_functions import _execute_r_script
from .etl import prepare_relative_abundance_data


def create_read_count_plots(melted_data, pct_data, output_file="read_count_tracking.html"):
    """Create Bokeh plots using iqplot for box plots"""
    
    sample_data = pct_data[pct_data['guess_type'] == 'sample']
    control_data = pct_data[pct_data['guess_type'] == 'control']

    # Plot 1: Read count tracking (absolute numbers)
    p1 = iqplot.strip(
        data=melted_data,
        q='Reads',
        cats='Step',
        q_axis="y",
        width=600,
        height=400,
        marker_kwargs=dict(alpha=0.0),
        title="Read Count Tracking"
    )

    # Add jittered points
    source = ColumnDataSource(melted_data)
    p1.scatter(jitter('Step', 0.2, range=p1.x_range), 'Reads', source=source, 
             size=8, alpha=0.6, color='red')

    p1.xaxis.major_label_orientation = 1.2
    p1.yaxis.axis_label = "Reads"

    # Plot 2: Read count tracking (percentages) with grouped box plots
    p2 = iqplot.strip(
        data=pct_data,
        q='Reads', 
        cats='Step',
        q_axis="y",
        width=600,
        height=400,
        marker_kwargs=dict(alpha=0.0),
        title="Read Count Tracking (Percentage of Input)"
    )


    if len(sample_data) > 0:
        source_sample = ColumnDataSource(sample_data)
        p2.scatter(jitter('Step', 0.15, range=p2.x_range), 'Reads', source=source_sample,
                 size=8, alpha=0.6, color='blue', legend_label='sample')

    if len(control_data) > 0:
        source_control = ColumnDataSource(control_data)
        p2.scatter(jitter('Step', 0.15, range=p2.x_range), 'Reads', source=source_control,
                 size=8, alpha=0.6, color='red', legend_label='control')

    p2.xaxis.major_label_orientation = 1.2
    p2.yaxis.axis_label = "Fraction of Input Reads"
    p2.legend.location = "top_right"

    # Save plots
    save(column(p1, p2), filename=output_file)
    print(f"Plots saved as '{output_file}'")

def generate_quality_profile_plots(forward_reads, reverse_reads, dataset_name, output_dir, plot_count=3):
    """
    Generate quality profile plots using R/ggplot2
    Must be written in R to take advantage of dada functionality
    """
    os.makedirs(f"{output_dir}/imgs/{dataset_name}", exist_ok=True)

    fwd_file = f"{output_dir}/imgs/{dataset_name}/read_quality_fwd.png"
    rvr_file = f"{output_dir}/imgs/{dataset_name}/read_quality_rvr.png"
    
    forward_read_list = ', '.join([f'"{f}"' for f in forward_reads[:min(plot_count, len(forward_reads))]])
    reverse_read_list = ', '.join([f'"{f}"' for f in reverse_reads[:min(plot_count, len(reverse_reads))]])
    
    r_script = f"""
    library(dada2)
    library(ggplot2)
    library(Rcpp)

    forward_reads <- c({forward_read_list})
    reverse_reads <- c({reverse_read_list})
    
    # Check if files exist and have content
    for(file in c(forward_reads, reverse_reads)) {{
        if(!file.exists(file) || file.size(file) == 0) {{
            cat("Warning: File", file, "is missing or empty\\n")
        }}
    }}

    # Forward reads plot
tryCatch({{
    p1 <- plotQualityProfile(forward_reads) + 
      geom_vline(xintercept = 240, col='black', lty=2) +
      geom_hline(yintercept = 30, col='grey60', lty=3)
    ggsave("{fwd_file}", plot=p1, dpi=320)
}}, error = function(e) {{
        cat("Error plotting forward reads:", e$message, "\\n")
    }})

    # Reverse reads plot  
tryCatch({{
    p2 <- plotQualityProfile(reverse_reads) + 
      geom_vline(xintercept = 200, col='red', lty=2) +
      geom_hline(yintercept = 30, col='grey60', lty=3)
    ggsave("{rvr_file}", plot=p2, dpi=320)
      }}, error = function(e) {{
        cat("Error plotting reverse reads:", e$message, "\\n")
    }})
    """
    
    _execute_r_script(r_script)

def create_contamination_plot(contamination_df, output_file):
    """Create contamination analysis plot"""
    
    # Create plot
    p = figure(width=800, height=400, title="Contamination Analysis",
               x_range=list(contamination_df['Sample'].unique()))

    for sample_name, sample_df in contamination_df.groupby("Sample"):
      bottom = 0
      for index_id, sample_data in sample_df.iterrows():
        height = sample_data.relabund      
        color = 'red' if sample_data.decontam else 'black'

        label = f"Contaminant: {sample_data.decontam}"
        p.vbar(x=[sample_name], top=[bottom + height], bottom=[bottom],
                   width=0.8, color=color, alpha=0.8, legend_label=label)
        bottom += height
    
    p.xaxis.major_label_orientation = 1.5
    p.yaxis.axis_label = "Relative Abundance (%)"
    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    
    save(p, filename=output_file)
    print(f"Contamination plot saved as {output_file}")

def create_stackbar_plot(relative_long_df, taxonomic_level, output_file, colors=None, title_suffix="", filter_pattern=None, sample_subset=None, min_abundance=1, title=None):
    """
    Create a single stacked bar plot for a given taxonomic level
    
    Parameters:
    relative_long_df: DataFrame in long format with relative abundances and taxonomy
    taxonomic_level: Taxonomic level to plot (e.g., 'phylum', 'genus', 'species')
    output_file: Path for output HTML file
    title_suffix: Additional text for plot title
    filter_pattern: Optional pattern to filter taxa
    sample_subset: List of sample names to include in plot (if None, includes all samples)
    min_abundance: Minimum relative abundance (%) to include taxa individually (others lumped into "Other")
    """
    if colors is None:
        colors=Category20[20]

    # Filter by sample subset if specified
    if sample_subset is not None:
        plot_df = relative_long_df[relative_long_df['sample'].isin(sample_subset)].copy()
    else:
        plot_df = relative_long_df.copy()
    
    # Apply taxonomic filter if specified
    if filter_pattern:
        plot_df = plot_df[plot_df[taxonomic_level].str.contains(filter_pattern, case=False, na=False)]
        title_suffix = f" ({filter_pattern})" if not title_suffix else title_suffix
    
    # Aggregate by taxonomic level and sample
    agg_df = plot_df.groupby(['sample', taxonomic_level])['relabund'].sum().reset_index()
    
    if min_abundance > 0:
        taxon_totals = agg_df.groupby(taxonomic_level)['relabund'].sum()
        abundant_taxa = taxon_totals[taxon_totals >= min_abundance].index.tolist()
        
        # Create "Other" category for low-abundance taxa
        agg_df.loc[~agg_df[taxonomic_level].isin(abundant_taxa), taxonomic_level] = "Other"
        
        # Re-aggregate after lumping
        agg_df = agg_df.groupby(['sample', taxonomic_level])['relabund'].sum().reset_index()
    

    # Get unique samples and taxa (preserve order if sample_subset provided)
    if sample_subset is not None:
        samples = [s for s in sample_subset if s in agg_df['sample'].unique()]
    else:
        samples = agg_df['sample'].unique()
    
    taxa = agg_df[taxonomic_level].unique()
    taxa = [t for t in taxa if pd.notna(t) and str(t).strip() != '']  # Remove empty/NaN

    # Generate colors
    if len(taxa) <= 20:
        colors = colors[:len(taxa)] if len(taxa) > 3 else ['#006666', '#FFCC33', '#724419'][:len(taxa)]
    else:
        colors = sns.color_palette("husl", len(taxa)).as_hex()
    
    color_map = dict(zip(taxa, colors))
    if "Other" in color_map:
        color_map["Other"] = "#808080"  # Grey
    
    # Pivot to wide format: one row per sample, one column per taxon.
    wide_df = (
        agg_df.pivot_table(index='sample', columns=taxonomic_level,
                           values='relabund', fill_value=0)
        .reindex(samples)
        .reindex(columns=taxa, fill_value=0)
        .reset_index()
    )
    source = ColumnDataSource(wide_df)

    hover = HoverTool(tooltips=[
          ("Sample", "@sample"),
          (taxonomic_level.title(), "$name"),
          ("Abundance", "@$name{0.00}%")
      ])
    # Create figure
    p = figure(
        x_range=samples,
        width=max(800, len(samples) * 80),  # Dynamic width based on sample count
        height=800,
        title=f"Relative Abundance - {taxonomic_level.title()}{title_suffix}",
        toolbar_location="above",
        tools=[hover, "pan", "wheel_zoom", "box_zoom", "reset", "save"]

    )
    p.vbar_stack(
        taxa,
        x='sample',
        width=0.8,
        source=source,
        color=[color_map[t] for t in taxa],
        alpha=0.8,
        legend_label=[str(t)[:30] for t in taxa],
    )
    
    # Styling
    p.xaxis.major_label_orientation = 1.2
    p.yaxis.axis_label = "Relative Abundance (%)"
    p.legend.label_text_font_size = "8pt"
    p.legend.nrows = 20
    
    leg = p.legend[0]
    p.add_layout(leg,'right')
    p.legend.click_policy = "hide"


    
    # Save plot
    save(p, filename=output_file, title=title)
    
    return p

def replot_relative_abundance(counts_file_path, taxonomy_file_path, output_dir,
                              taxonomic_levels=['phylum', 'genus', 'species'],
                              filter_patterns=None, sample_subset=None,
                              colors=None, dataset_name=None):
    """
    Regenerate relative-abundance stackbar plots from persisted pipeline
    outputs without rerunning the pipeline.

    Reads the counts table (e.g. ASVs_counts.csv or the decontaminated
    clean counts file) and ASV_taxonomy.csv, rebuilds the long-format
    relative abundance dataframe, and rewrites the HTML plots.
    """
    relative_long_df = prepare_relative_abundance_data(
        counts_file_path, taxonomy_file_path
    )
    return create_relative_abundance_stackbars(
        relative_long_df, output_dir,
        taxonomic_levels=taxonomic_levels,
        filter_patterns=filter_patterns,
        sample_subset=sample_subset,
        colors=colors,
        dataset_name=dataset_name,
    )

def create_relative_abundance_stackbars(relative_long_df, output_dir, 
                                       taxonomic_levels=['phylum', 'genus', 'species'],
                                       filter_patterns=None, sample_subset=None, colors=None, dataset_name=None):
    """
    Create relative abundance stacked bar plots at different taxonomic levels
    
    Parameters:
    relative_long_df: df with relative abundance data at all taxonomic levels
    output_dir: Output directory for plots
    taxonomic_levels: List of taxonomic levels to plot
    filter_patterns: Dict of {level: pattern} to filter specific taxa
    """
    
    plots_created = []
    
    # Create plots for each taxonomic level
    for level in taxonomic_levels:
        if level in relative_long_df.columns:
            
            # Standard plot
            output_file = f"{output_dir}/imgs/relative_abundance_{level}.html"
            title = f"{dataset_name}, {level}"
            create_stackbar_plot(relative_long_df, level, output_file, sample_subset=sample_subset, colors=colors, title=title)
            plots_created.append(output_file)
            
            # Filtered plot if pattern specified
            if filter_patterns and level in filter_patterns:
                pattern = filter_patterns[level]
                filtered_output_file = f"{output_dir}/relative_abundance_{level}_{pattern.lower()}.html"
                title = f"{title}, {pattern.lower()}"
                create_stackbar_plot(relative_long_df, level, filtered_output_file, filter_pattern=pattern, sample_subset=sample_subset, colors=colors, title=title)
                plots_created.append(filtered_output_file)
    
    return plots_created

def create_multiple_abundance_plots(relative_long_df, output_dir):
    """Create standard set of relative abundance plots"""
    
    # Standard plots
    plots = create_relative_abundance_stackbars(
        relative_long_df, 
        output_dir,
        taxonomic_levels=['phylum', 'genus', 'species']
    )
    
    # Filtered plots for specific groups
    filtered_plots = create_relative_abundance_stackbars(
        relative_long_df, 
        output_dir,
        taxonomic_levels=['species', 'genus'],
        filter_patterns={'species': 'cyano', 'genus': 'cyano'}
    )
    
    all_plots = plots + filtered_plots
    print(f"Created {len(all_plots)} relative abundance plots")
    return all_plots
