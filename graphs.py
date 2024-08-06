import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase
import numpy as np
import matplotlib.ticker as ticker

#helper functions

def round_up_to_next_significant(number):
    magnitude = 10 ** math.floor(math.log10(number))
    rounded_number = math.ceil(number / magnitude) * magnitude
    return rounded_number

def check_accession(acc, accessions, strand, strands):
    n = accessions.index(acc)
    for i in range(1, 6):
        x = n - i
        y = n + i
        if x >= 0:
            sx = strands[x]
            accx = accessions[x]
        else:
            sx = -10
            accx = '0'
        if y < len(accessions):
            sy = strands[y]
            accy = accessions[y]
        else:
            sy = -10
            accy = '0'
        if 0 < int(acc) - int(accx) < 6 and sx == strand:
            return True
        if 0 < int(accy) - int(acc) < 6 and sy == strand:
            return True
    return False

def format_k(x, pos):
    return f'{float(x / 1000000)}'
#main fucntions

def plot_gene_distribution(df, gpa, output_dir):
    """
    Visualize the distribution of genes along contigs and save each plot.

    Parameters:
    - df: DataFrame containing gene distribution data including Contig Length,
          Gene Count, Gene Size, and Strand information.
    - output_dir: Directory to save the plots.
    """
  
    

    min_value = 3
    max_value = 7
    
    # Define the colors for each value
    value_color_mapping = {
        3: '#030764',
        4: '#4B0082',
        5: '#8a4fd1',
        6: '#cb7bff',
        7: '#ff7e7e'
    }

    # Create a custom colormap based on the value-color mapping
    cmap = ListedColormap([value_color_mapping[i] for i in range(min_value, max_value + 1)])

    # Plot gene presence/absence for each contig
    for contig, row in df.iterrows():
        fig, ax = plt.subplots(figsize=(10, 5))

        # Plot genes for the contig
        gene_starts = [int(x) for x in row['Gene Start']]
        gene_sizes = [int(x) for x in row['Gene Size']]
        gene_strands = row['Strand']
        gene_products = row['Gene Product']
        

        count = 0
        for start, size, strand, product in zip(gene_starts, gene_sizes, gene_strands, gene_products):
            count += 1
            num = sum(gpa.loc[product])
            color = cmap(num - min_value)  # Adjusted the mapping

            marker = 'o'
            if 'nitrate reductase' in product.lower():
                marker = 'x'

            ax.plot(start, count, color=color, marker=marker, markersize=1 + size / 500)

        # Set y-axis ticks and labels
        ax.set_yticks([x for x in range(1, count + 1)])
        ax.set_xlabel('Position in Mbp')
        ax.set_ylabel('Genes')
        #ax.set_title(f"Distribution of exclusive genes along contig {contig}")

        # Set x-axis ticks and labels
        contig_length = round_up_to_next_significant(int(row['Contig Length']))
        ax.set_xticks(range(0, contig_length, 100000))
        ax.set_xticklabels([format_k(x, None) if i % 5 == 0 else '' for i, x in enumerate(ax.get_xticks())])
        ax.set_yticks(range(0, count, 10))
        ax.set_yticklabels(range(0, count, 10))
        
        ax.grid(True)
        ax.grid(which='minor', axis='x', linestyle=':', linewidth='0.5', color='gray')
        ax.minorticks_on()

        # Major gridlines: every 5th column
        ax.set_xticks(range(0, contig_length, 500000), minor=False)
        ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
        # Create colorbar legend
        cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
        cb = ColorbarBase(plt.gca(), cmap=cmap, orientation='vertical')
        cb.set_ticklabels([str(i) for i in range(3,8)])
        cb.set_label('Number of Genomes sharing gene')

        # Add a legend for the marker sizes
        marker_sizes = [1000, 2000, 3000, 4000 ]
        marker_labels = ['1000bp', '2000bp', '3000bp', '4000bp']
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=1 + size / 500) for size in marker_sizes]
        ax.legend(handles, marker_labels, title="Gene Size", loc='upper left', bbox_to_anchor=(0.01, 1))

        # Save the plot
        output_filename = f"{output_dir}/{contig}_gene_distribution.png"
        plt.savefig(output_filename)

        # Close the plot to release memory
        plt.close(fig)

def plot_gene_clusters(df,gpa, output_dir):
    
    min_value = 3
    max_value = 7
    
    # Define the colors for each value
    value_color_mapping = {
        3: '#030764',
        4: '#4B0082',
        5: '#8a4fd1',
        6: '#cb7bff',
        7: '#ff7e7e'
    }

    # Create a custom colormap based on the value-color mapping
    cmap = ListedColormap([value_color_mapping[i] for i in range(min_value, max_value + 1)])
    selected_rows = []
    # Plot gene presence/absence for each contig
    for contig, row in df.iterrows():
        fig, ax = plt.subplots(figsize=(10, 5))

        # Plot genes for the contig
        gene_starts = [int(x) for x in row['Gene Start']]
        gene_sizes = [int(x) for x in row['Gene Size']]
        gene_strands = row['Strand']
        gene_products = row['Gene Product']
        gene_accessions=row['Accession']

        filtered_starts = []
        filtered_sizes = []
        filtered_strands = []
        filtered_products = []
        filtered_accessions = []

        count = 0
        for start, size, strand, product, acc in zip(gene_starts, gene_sizes, gene_strands, gene_products, gene_accessions):
            count += 1
            num = sum(gpa.loc[product])
            color = cmap(num - min_value)  

            marker = 'o'
            if 'nitrate reductase' in product.lower():
                marker = 'x'
            if check_accession(acc, gene_accessions, strand, gene_strands): 
                ax.plot(start, count, color=color, marker=marker, markersize=1 + size / 500)
                filtered_starts.append(start)
                filtered_sizes.append(size)
                filtered_strands.append(strand)
                filtered_products.append(product)
                filtered_accessions.append(acc)

            filtered_gene_data = {
                'Gene Start': filtered_starts,
                'Gene Size': filtered_sizes,
                'Strand': filtered_strands,
                'Gene Product': filtered_products,
                'Accession': filtered_accessions
            }
        
        selected_df = pd.DataFrame(filtered_gene_data)
        selected_df.to_csv(f'{output_dir}/{contig}_gene_clusters.csv')
        # Set y-axis ticks and labels
        ax.set_yticks([x for x in range(1, count + 1)])
        ax.set_xlabel('Contig')
        ax.set_title(f"Exclusive gene clusters along contig {contig}")

        # Set x-axis ticks and labels
        contig_length = round_up_to_next_significant(int(row['Contig Length']))
        ax.set_xticks(range(0, contig_length, int(contig_length / 10)))
        ax.set_xticklabels(range(0, contig_length, int(contig_length / 10)))
        ax.set_yticks(range(0, count, 5))
        ax.set_yticklabels(range(0, count, 5))

        # Create colorbar legend
        cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
        cb = ColorbarBase(plt.gca(), cmap=cmap, orientation='vertical')
        cb.set_ticklabels([str(i) for i in range(3,8)])
        cb.set_label('Number of Genomes sharing gene')

        # Save the plot
        output_filename = f"{output_dir}/{contig}_gene_clusters.png"
        plt.savefig(output_filename)

        # Close the plot to release memory
        plt.close(fig)

def get_pa_matrix(gpa): 
    import rpy2
    from rpy2.robjects import pandas2ri, r
    from rpy2.robjects.packages import importr
    import pandas as pd
    
    # Activate pandas2ri conversion
    pandas2ri.activate()
    
    # Load the required R packages
    if not rpy2.robjects.packages.isinstalled("heatmaply"):
        utils = importr('utils')
        utils.install_packages("heatmaply")
    
    if not rpy2.robjects.packages.isinstalled("htmlwidgets"):
        utils.install_packages("htmlwidgets")
        
    heatmaply = importr('heatmaply')
    htmlwidgets = importr('htmlwidgets')
    
    # Read the CSV file into a pandas DataFrame
    
    # Convert pandas DataFrame to R DataFrame
    r_dataframe = pandas2ri.py2rpy(gpa)
    
    # Convert R DataFrame to R matrix
    r_matrix = r['as.matrix'](r_dataframe)
    # Define R code for generating the heatmap
    print(gpa.dtypes)
    r_code = f'''
    r_matrix <- {r_matrix.r_repr()}
    print(r_matrix)
    heatmap <- heatmaply(r_matrix, labRow = FALSE, scale_fill_binary = TRUE)
    htmlwidgets::saveWidget(heatmap, "heatmap_output.html")
    '''.format(r_matrix.r_repr())
    
    # Execute R code
    r(r_code)
