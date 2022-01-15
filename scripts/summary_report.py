#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as offline
import argparse
from pathlib import Path
import numpy as np


#need to pip install -U kaleido

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_dir', metavar='spear_vcfs/merged.spear.vcf', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('output_filename', metavar='spear_vcfs/merged.spear.vcf', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS
    parser.add_argument('baseline', metavar='Omicron', type=str,
        help='lineage for baseline') #ADD A DEFAULT FOR THIS
    args = parser.parse_args()

    scores_summary = pd.read_csv(f'{args.input_dir}/spear_score_summary.tsv', sep = '\t')
    annotation_summary = pd.read_csv(f'{args.input_dir}/spear_annotation_summary.tsv', sep = '\t')
    annotation_summary["compound_nt_var"] = annotation_summary["REF"] + annotation_summary["POS"].astype("str") + annotation_summary["ALT"]
    baseline_scores = pd.read_csv(f'~/Desktop/multi_omi_new_summarise/spear_score_summary.tsv', sep = '\t')

    variants_counts_table = annotation_summary.groupby(["REF","POS", "ALT","compound_nt_var"])[["compound_nt_var"]].count()
    variants_counts_table.columns = ["count"]
    variants_counts_table = variants_counts_table.sort_values("count", axis = 0 , ascending = False)
    variants_counts_table = variants_counts_table.reset_index()

    variants_table = go.Figure(data=[go.Table(
        header=dict(values= ["POS","REF", "ALT", "count"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[variants_counts_table["REF"], variants_counts_table["POS"], variants_counts_table["ALT"], variants_counts_table["count"]],
                   fill_color='lavender',
                   align='left'))
    ])
    variants_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)'})
    variants_table.write_html("variants_table.html")   
    variants_table_plt = offline.plot(variants_table,output_type='div', include_plotlyjs = True, config = {'displaylogo': False}) #including plotly js with this plot and not with any future ones to keep html size down 

    residues_counts_table = annotation_summary.groupby("residues")[["residues"]].count()
    residues_counts_table.columns = ["count"]
    residues_counts_table = residues_counts_table.sort_values("count", axis = 0, ascending = False)
    residues_counts_table = residues_counts_table.reset_index()
    
    residues_table = go.Figure(data=[go.Table(
        header=dict(values=list(residues_counts_table.columns),
                    fill_color='paleturquoise'),
                    #align='left'),
        cells=dict(values=[residues_counts_table["residues"], residues_counts_table["count"]],
                   fill_color='lavender'))
                   #align='left'))
    ])  
    residues_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)'})
    residues_table.write_html("residue_table.html")
    residues_table_plt = offline.plot(residues_table,output_type='div', include_plotlyjs = False , config = {'displaylogo': False})
    
    graph_scores = baseline_scores.columns.tolist()
    non_graph_scores = ['sample_id', 'total_variants', 'total_residue_variants','consequence_type_variants', 'region_residues', 'domain_residues','ACE2_contact_counts', 'ACE2_contact_score', 'trimer_contact_counts','trimer_contact_score', 'barns_class_variants']
    graph_scores = [score for score in graph_scores if score not in non_graph_scores]
    baseline_scores_array = baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, graph_scores].fillna(0).values.squeeze()
    sample_graph_scores = scores_summary[graph_scores]
    sample_graph_ids = scores_summary["sample_id"].values.tolist()
    sample_baseline_relative_scores = sample_graph_scores - baseline_scores_array

    cols_to_display = sample_baseline_relative_scores.columns[sample_baseline_relative_scores.isnull().all() == False].values.tolist()
    horizontal_bar_plot = go.Figure()
    for column in cols_to_display:
        horizontal_bar_plot.add_trace(
            go.Bar(
                x = sample_baseline_relative_scores[column],
                y = sample_graph_ids,
                text = sample_graph_ids,
                orientation = 'h',
                name = column,
                textposition = "inside"
            )
        )
    buttons = []
    count = 0
    for col in cols_to_display:
        trace_list = [False] * len(cols_to_display)
        trace_list[count] = True
        buttons.append(
            dict(
                label = col,
                method = 'update',
                args = [
                    {'visible': trace_list},
                    {'title': col,'showlegend':True}])
            )
        count += 1
    horizontal_bar_plot.update_layout(
    updatemenus=[go.layout.Updatemenu(
        active=0,
        buttons = buttons
        )
    ])

    horizontal_bar_plot.update_layout({"paper_bgcolor":'rgba(0,0,0,0)'})
    horizontal_bar_plot_html = offline.plot(horizontal_bar_plot, output_type = 'div', include_plotlyjs=False ,config = {'displaylogo': False, 'displayModeBar': False})
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <nav class="navbar navbar-expand-lg py-3 navbar-dark bg-dark shadow-sm">
                <div class="container">
                    <a href="#" class="navbar-brand">
                    <!-- Logo Image -->
                    <img src="SPEAR_smallest.png" class="img-fluid" alt="Responsive image" width="100" class="d-inline-block align-middle mr-2">
                    <!-- Logo Text -->
                    <span class="text-uppercase font-weight-bold">SPEAR Summary Report</span>
                    </a>
                </div>
            </nav>
            <div class="container">
                <div class="row">
                    <h2>Section 1: Variant and residue counts</h2>
                </div>
                <div class="row">
                    <div class="col-6">
                        ''' + variants_table_plt + '''
                        <div id="plot1" class="responsive-plot"></div>
                    </div>
                    <div class="col-6">
                        ''' + residues_table_plt + '''
                        <div id="plot1" class="responsive-plot"></div>
                    </div>
                </div>
                <div class = "row">
                    <div class="col-md">
                        <p>Summary of nucleotide variants and residue changes across all samples</p>
                    </div>
                </div>
                <div class = "row">
                    <div class = "col-12">
                        <h2>Scoring Summary</h2>
                    </div>
                </div>
                <div class = "row">
                    <div class = "col-12">
                        ''' + horizontal_bar_plot_html + '''
                    </div>
                </div>
                <div class = "row">
                    <div class = "col-12">
                        <p>Plot of scores, use dropdown menu to view individual scores for all samples.</p>
                    </div>
                </div>
            </div> 
        <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
        <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
        </body>
    </html>'''

    f = open('report.html','w')
    f.write(html_string)
    f.close()

    

# this one lets you colour cells by variable 
#     import plotly.graph_objects as go
# from plotly.colors import n_colors
# import numpy as np
# np.random.seed(1)

# colors = n_colors('rgb(255, 200, 200)', 'rgb(200, 0, 0)', 9, colortype='rgb')
# a = np.random.randint(low=0, high=9, size=10)
# b = np.random.randint(low=0, high=9, size=10)
# c = np.random.randint(low=0, high=9, size=10)

# fig = go.Figure(data=[go.Table(
#   header=dict(
#     values=['<b>Column A</b>', '<b>Column B</b>', '<b>Column C</b>'],
#     line_color='white', fill_color='white',
#     align='center',font=dict(color='black', size=12)
#   ),
#   cells=dict(
#     values=[a, b, c],
#     line_color=[np.array(colors)[a],np.array(colors)[b], np.array(colors)[c]],
#     fill_color=[np.array(colors)[a],np.array(colors)[b], np.array(colors)[c]],
#     align='center', font=dict(color='white', size=11)
#     ))
# ])

# fig.show()


#sub plots

# from plotly.subplots import make_subplots
# import plotly.graph_objects as go

# fig = make_subplots(rows=1, cols=2)

# fig.add_trace(
#     go.Scatter(x=[1, 2, 3], y=[4, 5, 6]),
#     row=1, col=1
# )

# fig.add_trace(
#     go.Scatter(x=[20, 30, 40], y=[50, 60, 70]),
#     row=1, col=2
# )

# fig.update_layout(height=600, width=800, title_text="Side By Side Subplots")
# fig.show()

if __name__ == "__main__":
    main()