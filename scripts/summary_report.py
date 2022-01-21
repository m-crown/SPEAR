#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as offline
import argparse
from pathlib import Path
import numpy as np
from shutil import copyfile
from numpy.random import seed
from numpy.random import uniform



#need to pip install -U kaleido

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('score_summary', metavar='spear_score_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS
    parser.add_argument('annotation_summary', metavar='spear_annotation_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('baseline_scores', metavar='baseline.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('images_dir', metavar='spear_images/', type=str,
        help='Directory for spear images to be copied from') 
    parser.add_argument('output_dir', metavar='report/', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS
    parser.add_argument('baseline', metavar='Omicron', type=str,
        help='lineage for baseline') #ADD A DEFAULT FOR THIS
    args = parser.parse_args()

    Path(f'{args.output_dir}/images').mkdir(parents=True, exist_ok=True)

    scores_summary = pd.read_csv(f'{args.score_summary}', sep = '\t')
    #print(scores_summary)
    annotation_summary = pd.read_csv(f'{args.annotation_summary}', sep = '\t')
    annotation_summary["compound_nt_var"] = annotation_summary["REF"] + annotation_summary["POS"].astype("str") + annotation_summary["ALT"]
    #print(annotation_summary)
    baseline_scores = pd.read_csv(f'{args.baseline_scores}', sep = '\t')
    #print(baseline_scores)

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
    residues_table_plt = offline.plot(residues_table,output_type='div', include_plotlyjs = False , config = {'displaylogo': False})
    
    graph_scores = baseline_scores.columns.tolist()
    non_graph_scores = ['sample_id', 'total_variants', 'total_residue_variants','consequence_type_variants', 'region_residues', 'domain_residues','ACE2_contact_counts', 'ACE2_contact_score', 'trimer_contact_counts','trimer_contact_score', 'barns_class_variants','bloom_ACE2_max', 'bloom_ACE2_min', 'VDS_max', 'VDS_min', 'serum_escape_max', 'serum_escape_min', 'mAb_escape_all_classes_max', 'mAb_escape_all_classes_min', 'mAb_escape_class_1_max', 'mAb_escape_class_1_min', 'mAb_escape_class_2_max', 'mAb_escape_class_2_min', 'mAb_escape_class_3_max', 'mAb_escape_class_3_min', 'mAb_escape_class_4_max', 'mAb_escape_class_4_min', 'BEC_RES_max', 'BEC_RES_min', 'BEC_EF_sample']
    graph_scores = [score for score in graph_scores if score not in non_graph_scores]
    baseline_scores_array = baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, graph_scores].fillna(0).values.squeeze()
    

    sample_graph_scores = scores_summary[graph_scores]
    sample_scores = sample_graph_scores.dropna(axis=1, how = "all") #remove empty cols from table to be displayed (do this later for the graph table to allow subtraction of baseline array)
    sample_graph_ids = scores_summary["sample_id"].values.tolist()
    #print("HELLO" , baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, "sample_id"])
    sample_table_ids = baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, "sample_id"].values.tolist() + scores_summary["sample_id"].values.tolist()
    #print(sample_table_ids)
    sample_baseline_relative_scores = sample_graph_scores - baseline_scores_array
    sample_scores = pd.concat([baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, graph_scores].fillna(0), sample_scores])
    #print(sample_scores)
    cols = sample_baseline_relative_scores.columns[sample_baseline_relative_scores.isnull().all() == False].values.tolist()
    cols_to_display = sample_baseline_relative_scores.columns[sample_baseline_relative_scores.isnull().all() == False].str.replace("_", " ").values.tolist()
    sample_scores_baseline_relative = np.greater(sample_baseline_relative_scores[cols].to_numpy(), 0) #all others should be coloured if greater than baseline
    samples_scores_baseline_colours = np.where(sample_scores_baseline_relative == True, "rgb(250,180,174)", "rgb(179,205,227)")
    scores_table = go.Figure(data=[go.Table(
        header=dict(values= ["Sample ID"] + cols_to_display,
                    fill_color='paleturquoise',
                    align='center'),
        cells=dict(values=[sample_table_ids] + [sample_scores[x].round(2).fillna("").astype("str") for x in cols],
                   fill_color=[['lavender'] * len(sample_graph_ids)] + samples_scores_baseline_colours.T.tolist(),
                   align='center'))
    ])  
    scores_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)'})
    scores_table_plt = offline.plot(scores_table, output_type='div', include_plotlyjs = False , config = {'displaylogo': False})

    horizontal_bar_plot = go.Figure()
    colors = [
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#d62728',  # brick red
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
    ]

    for column in cols:
        horizontal_bar_plot.add_trace(
            go.Bar(
                x = sample_baseline_relative_scores[column],
                y = sample_graph_ids,
                text = sample_graph_ids,
                orientation = 'h',
                name = column,
                textposition = "inside",
                marker=dict(
                    color=colors[cols.index(column)],
                    line=dict(color='rgb(248, 248, 249)', width=1)
            )
                
            )
        )
    buttons = []
    trace_list = [True] * len(cols_to_display)
    buttons.append(
            dict(
                label = "All Scores",
                method = 'update',
                args = [
                    {'visible': trace_list},
                    {'title': "All",'showlegend':True}])
            )
    count = 0
    cols_to_display = ["All"] + cols_to_display
    for col in cols_to_display:
        if col == "All":
            continue
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
    horizontal_bar_plot_html = offline.plot(horizontal_bar_plot, output_type = 'div', include_plotlyjs=False ,config = {'displaylogo': False})
    
    #use the pkl file of all ORFs to make a heatmap, where you merge the respos from anno df with 
    respos_df = pd.read_csv(f'data/product_mapping.csv')
    #print(respos_df)
    #print(annotation_summary)
    fig = px.bar(respos_df, x="residue", y="organism", color="product", title="Long-Form Input", orientation = "h")

    #for insertions these are annotated as a change on the last ref res 
    annotation_summary['refres'] = annotation_summary["residues"].str.extract('([A-Z])-*[0-9]+-*[a-zA-Z]+')
    annotation_summary['respos'] = annotation_summary["residues"].str.extract('[A-Z]-*([0-9]+)-*[a-zA-Z]+')
    annotation_summary['altres'] = annotation_summary["residues"].str.extract('[A-Z]-*[0-9]+-*([a-zA-Z]+)')
    #print(annotation_summary.loc[annotation_summary["refres"].isna()])
    #print(annotation_summary.loc[annotation_summary["residues"].fillna("").str.contains("EPE", regex = True)])
    
    occurences = annotation_summary.groupby(['description' , 'refres','respos', "residues"]).size().reset_index().rename(columns={0:'count'})
    #y = {}
    #occurences['y'] = occurences['description'].map(y) 
    occurences["respos"] = occurences["respos"].astype(int)
    
    orf1ab = ['leader', 'nsp2', 'nsp3', 'nsp4', '3C-like proteinase', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'nsp11', 'RNA-dependent RNA polymerase', 'helicase', "3'-to-5' exonuclease", 'endoRNAse', "2'-O-ribose methyltransferase"]
    structural = ['surface glycoprotein','nucleocapsid phosphoprotein', 'envelope protein', 'membrane glycoprotein']
    accessory = ['ORF3a protein', 'ORF6 protein', 'ORF7a protein', 'ORF7b', 'ORF8 protein', 'ORF10 protein']
    products = orf1ab + structural + accessory
    product_colours = {
        'leader': "lightsteelblue", 
        'nsp2': "lightskyblue", 
        'nsp3':"lightsteelblue", 
        'nsp4':"lightskyblue", 
        '3C-like proteinase':"lightsteelblue", 
        'nsp6':"lightskyblue", 
        'nsp7':"lightsteelblue", 
        'nsp8':"lightskyblue", 
        'nsp9':"lightsteelblue", 
        'nsp10':"lightskyblue", 
        'nsp11':"lightsteelblue", 
        'RNA-dependent RNA polymerase':"lightskyblue", 
        'helicase':"lightsteelblue", 
        "3'-to-5' exonuclease":"lightskyblue", 
        'endoRNAse':"lightsteelblue", 
        "2'-O-ribose methyltransferase":"lightskyblue",
        'surface glycoprotein': "red",
        'nucleocapsid phosphoprotein': "green", 
        'envelope protein': "lightblue", 
        'membrane glycoprotein': "darkblue",
        'ORF3a protein': "orange",
        'ORF6 protein': "gold", 
        'ORF7a protein': "brown", 
        'ORF7b': "purple", 
        'ORF8 protein': "cyan",
        'ORF10 protein': "magenta"}
    merged_occurences = pd.merge(occurences, respos_df, left_on = ["description", "respos"], right_on = ["product", "residue-position"], how = "left")
    prods_in_samples = [x for x in products if x in merged_occurences["product"].unique()]

    axis_products = {val:prods_in_samples.index(val) for val in prods_in_samples}
    merged_occurences['axis_pos'] = merged_occurences['product'].map(axis_products)
    merged_orf_occurences = merged_occurences.loc[merged_occurences["residue"].isna() == False]
    merged_orf_occurences["product"] = merged_orf_occurences["product"].astype("category")
    merged_orf_occurences["product"] = merged_orf_occurences["product"].cat.set_categories(products)
    merged_orf_occurences = merged_orf_occurences.sort_values(["product", "respos"])
    shapes_to_plot = merged_occurences[["product", "axis_pos"]].groupby(["product", "axis_pos"]).size().reset_index()
    shapes_to_plot["min_pos"] = shapes_to_plot["product"].apply(lambda x: respos_df.loc[respos_df['product'] == x, "position"].min())
    shapes_to_plot["max_pos"] = shapes_to_plot["product"].apply(lambda x: respos_df.loc[respos_df['product'] == x, "position"].max())
    shapes_to_plot["axis_pos"] = 0
    seed(42069)
    positive = True
    for prod in prods_in_samples:
        if prod in orf1ab:
            if positive:
                merged_orf_occurences.loc[merged_orf_occurences["product"] == prod, "axis_pos"] = np.linspace(start = 5, stop = 1, num = len(merged_orf_occurences.loc[merged_orf_occurences["product"] == prod]))
                merged_orf_occurences.loc[merged_orf_occurences["product"] == prod, "textpos"] = "top center"
                positive = False
            else:
                merged_orf_occurences.loc[merged_orf_occurences["product"] == prod, "axis_pos"] = np.linspace(start = -5, stop = -1, num = len(merged_orf_occurences.loc[merged_orf_occurences["product"] == prod]))
                merged_orf_occurences.loc[merged_orf_occurences["product"] == prod, "textpos"] = "bottom center"
                positive = True
    #rest are direct
    merged_orf_occurences.loc[merged_orf_occurences["product"] == "surface glycoprotein", "axis_pos"] = np.linspace(start = 10.5, stop = 5.5, num = len(merged_orf_occurences.loc[merged_orf_occurences["product"] == "surface glycoprotein"]))
    merged_orf_occurences.loc[merged_orf_occurences["product"] == "surface glycoprotein", "textpos"] = "top center"
    merged_orf_occurences.loc[merged_orf_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"]), "axis_pos"] = np.linspace(start = 10.5, stop = 5.5, num = len(merged_orf_occurences.loc[merged_orf_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"])]))
    merged_orf_occurences.loc[merged_orf_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"]), "textpos"] = "top center"
    merged_orf_occurences.loc[merged_orf_occurences["product"].isin(accessory), "axis_pos"] = np.linspace(start = -5, stop = -10, num = len(merged_orf_occurences.loc[merged_orf_occurences["product"].isin(accessory)]))
    merged_orf_occurences.loc[merged_orf_occurences["product"].isin(accessory), "textpos"] = "bottom center"
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=merged_orf_occurences["position"], y=merged_orf_occurences["axis_pos"], text = merged_orf_occurences["residues"],mode = "markers+text", textposition= merged_orf_occurences["textpos"]))
    fig.update_traces(
        marker={
            "size" : 12,
            "line" : {"width" : 2 , "color": 'DarkSlateGrey'}},
        selector=dict(mode='markers+text'))
    fig.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',xaxis=dict(showgrid=False, showline = False, zeroline = False),yaxis=dict(showgrid=False, showline = False, zeroline = False)) #title_text='WHAT SHOULD THIS TITLE BE ? ', title_x=0.5, 
    def add_orf_shape(x_min_pos, x_max_pos, product, orf1ab, structural, accessory,product_colours):
        shapes = []
        if product in orf1ab:
            y_start = 0
        elif product in structural:
            y_start = 4
        elif product in accessory:
            y_start = -4
        else:
            y_start = 10 #shouldnt ever be set but just in case it makes it obvious
        x0 = x_min_pos
        y0 = y_start - 0.5
        x1 = x_max_pos
        y1 = y_start + 0.5
        h_x = min((x1-x0)*0.2,100) #(x1-x0) * 0.2
        h_y = (y1-y0) *0.2
        rounded_bottom_left = f' M {x0+h_x}, {y0} Q {x0}, {y0} {x0}, {y0+h_y}'
        rounded_top_left = f' L {x0}, {y1-h_y} Q {x0}, {y1} {x0+h_x}, {y1}'
        rounded_top_right = f' L {x1-h_x}, {y1} Q {x1}, {y1} {x1}, {y1-h_y}'
        rounded_bottom_right = f' L {x1}, {y0+h_y} Q {x1}, {y0} {x1-h_x}, {y0}Z'
        path = rounded_bottom_left + rounded_top_left+\
                rounded_top_right+rounded_bottom_right
        shape = dict(
            path=path,
            fillcolor= product_colours[product],
            layer='above', 
            line=dict(color='rgba(128, 0, 128, 0.7)', width=0.5))
        return shape
    #draw orf boxes from the respos_df using min and max of each group of groubpy product
    orf_boxes = shapes_to_plot.apply(lambda x: add_orf_shape(x.min_pos, x.max_pos, x["product"], orf1ab, structural, accessory, product_colours), axis =1).values.tolist()

    def add_mutation_line(x_pos, y_pos, prod, orf1ab, structural, accessory):
        if prod in orf1ab:
            y = 0
        elif prod in structural:
            y = 4
        elif prod in accessory:
            y = -4
        else:
            y = 10 #shouldnt ever be set but just in case it makes it obvious
        shape = dict(
            type="line",
            x0=x_pos, y0=y, x1=x_pos, y1=y_pos,
            line=dict(
                color="grey",
                width=2,
                ),
            layer='below'
            )
        return shape
    
    lines = merged_orf_occurences.apply(lambda x: add_mutation_line(x.position, x.axis_pos, x["product"], orf1ab, structural, accessory),axis = 1).values.tolist()
    combo = orf_boxes + lines
    fig.update_layout(shapes = combo)
    # Add dropdown
    fig.update_layout(
        updatemenus=[
            dict(
                buttons=list([
                    dict(
                        args=["mode", "markers+text"],
                        label="Show text",
                        method="restyle"),
                    dict(
                        args=["mode", "markers"],
                        label="Hide text",
                        method="restyle"
                    ),
                    dict(
                        args=["mode", "text"],
                        label="Text only",
                        method="restyle"
                    ),
                ]),
                active = 0,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.08,
                yanchor="top"
            ),
            dict(
                buttons=list([
                    dict(
                        args=[{"shapes" : combo}],
                        label="Show shapes",
                        method="relayout"),
                    dict(
                        args=[{"shapes" : orf_boxes}],
                        label="Hide shapes",
                        method="relayout"),
                ]),
                active = 0,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.2,
                xanchor="left",
                y=1.08,
                yanchor="top"
            ),
        ]
    )

    #fig.show()
    mutation_plot_html = offline.plot(fig, output_type = 'div', include_plotlyjs=False ,config = {'displaylogo': False})

    #going to make a line set for each sample: ORF positions and points can be dynamically controlled.
    #can have a line set for all and set data according to this. 
    print(merged_orf_occurences.loc[merged_orf_occurences["respos"] == 1033])
    print(merged_occurences[["respos", "residue"]].groupby("respos").count())
    #can do a cluttered plot at top (or near top of page for all samples) that just has a warning that for lots of samples it is very cluttered and to refer to the table of sample plot links. Or when single sample have a report with everything at top, when multi have a table and lots of plots. 
    

    #fig.write_image(f'{args.experiment}.png')
    #fig.update_traces(marker=dict(size=12,
    #                          line=dict(width=2,
    #                                    color='DarkSlateGrey')),
    #              selector=dict(mode='markers'))
    #fig.write_html(f'{args.experiment}.html')

    #fig.write_html("residues.html")
    #ORF dataframe with the residue change info from the annotation summary 
    #for each sample and the HDF5 file domain info to get a heatmap that on hover 
    # shows details about a residue
    
    #maybe a blobbogram of median min max and a baseline line? could have a dropdown on the horizontal bar plot to change the data vis style. 
    
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
                    <img src="images/SPEAR_smallest.png" class="img-fluid" alt="Responsive image" width="50" class="d-inline-block align-left mr-2">
                    <!-- Logo Text -->
                    <span class="text-uppercase font-weight-bold">SPEAR Summary Report</span>
                    </a>
                </div>
            </nav>
            <div class="container-fluid">
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
                    <div class="col-md">
                        ''' + scores_table_plt + '''
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
                <div class = "row">
                    <div class = "col-12 h-100">
                        ''' + mutation_plot_html + '''
                    </div>
                </div>
            </div> 
        <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
        <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
        </body>
    </html>'''
    
    copyfile(f'{args.images_dir}/SPEAR_smallest.png', f'{args.output_dir}/images/SPEAR_smallest.png')
    f = open(f'{args.output_dir}/report.html','w')
    f.write(html_string)
    f.close()

if __name__ == "__main__":
    main()