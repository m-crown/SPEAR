#!/usr/bin/env python3

from re import sub
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as offline
import plotly.io as io
import argparse
from pathlib import Path
import numpy as np
from shutil import copyfile
from numpy.random import seed

from rich.console import Console
from rich.table import Table
from rich.progress import track

#need to pip install -U kaleido

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
        line=dict(color='black', width=0.5))
    return shape

def add_orf_label(x_min_pos , x_max_pos, product, orf1ab, structural, accessory):
    x_mid_pos = x_max_pos - ((x_max_pos - x_min_pos)/2)
    product_mapping = {
        "leader protein" : "NSP1", 
        "nsp2" : "NSP2", 
        "nsp3" : "NSP3" , 
        "nsp4" : "NSP4", 
        "3C-like proteinase" : "NSP5", 
        "nsp6": "NSP6", 
        "nsp7" : "NSP7", 
        "nsp8" : "NSP8", 
        "nsp9" : "NSP9", 
        "nsp10" : "NSP10", 
        "nsp11" : "NSP11", 
        "RNA-dependent RNA polymerase" : "NSP12", 
        "helicase" : "NSP13", 
        "3'-to-5' exonuclease" : "NSP14",
        "endoRNAse" : "NSP15",
        "2'-O-ribose methyltransferase" : "NSP16", 
        'surface glycoprotein': "S",
        'nucleocapsid phosphoprotein': "N", 
        'envelope protein': "E", 
        'membrane glycoprotein': "M",
        'ORF3a protein': "3a",
        'ORF6 protein': "6", 
        'ORF7a protein': "7a", 
        'ORF7b': "7b", 
        'ORF8 protein': "8",
        'ORF10 protein': "10"}
    if product != "nsp11":
        prod = product_mapping[product]
    else:
        prod = "" #trying to fit this on 14AA box is impossible. 
    if product in orf1ab:
        y_pos = 0
    elif product in structural:
        y_pos = 4
    elif product in accessory:
        y_pos = -4
    else:
        y_pos = 10 #shouldnt ever be set but just in case it makes it obvious
    
    annotation = go.layout.Annotation(x = x_mid_pos, y = y_pos, text = prod, font = {'size': 10}, showarrow = False, align = "center")
    return annotation


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


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--product_plots', default=False, action='store_true',
        help = "Output per sample ORF plots.")
    parser.add_argument('score_summary', metavar='spear_score_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS
    parser.add_argument('annotation_summary', metavar='spear_annotation_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('baseline_scores', metavar='baseline.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS 
    parser.add_argument('input_samples', metavar='1000', type=str,
        help='Number of samples input into SPEAR pipeline') #ADD A DEFAULT FOR THIS 
    parser.add_argument('qc_samples', metavar='700', type=str,
        help='Number of samples passing QC and entering SPEAR pipeline') #ADD A DEFAULT FOR THIS 
    parser.add_argument('images_dir', metavar='spear_images/', type=str,
        help='Directory for spear images to be copied from') 
    parser.add_argument('scripts_dir', metavar='$CONDA_PREFIX/bin', type=str,
        help='Directory for spear scripts')
    parser.add_argument('data_dir', metavar='$CONDA_PREFIX/data', type=str,
        help='Directory for spear data')        
    parser.add_argument('output_dir', metavar='report/', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF') #ADD A DEFAULT FOR THIS
    parser.add_argument('baseline', metavar='Omicron', type=str,
        help='lineage for baseline') #ADD A DEFAULT FOR THIS
    args = parser.parse_args()

    #eventually make this an argument:
    cli_out = True
    #needs to take as input the conda bin location where plotly js is stored!
    #NEED TO GET THE MIN AND MAX OF EACH SCORE SCALE IN A FILE TO READ IN AND USE TO SET COLOR SCALES - OR JUST WORK IT OUT AND HARD CODE IT ? 
    #need to make a function for each of the manipulations of dataframes to streamline it and allow things like changing baseline? ? 
    seed(42069)
    
    Path(f'{args.output_dir}/images').mkdir(parents=True, exist_ok=True)
    Path(f'{args.output_dir}/plots/plotly').mkdir(parents=True, exist_ok=True)
    copyfile(f'{args.scripts_dir}/plotly-2.8.3.min.js', f'{args.output_dir}/plots/plotly/plotly-2.8.3.min.js')
    Path(f'{args.output_dir}/plots/product_plots').mkdir(parents=True, exist_ok=True)

    scores_summary = pd.read_csv(f'{args.score_summary}', sep = '\t')
    annotation_summary = pd.read_csv(f'{args.annotation_summary}', sep = '\t')
    annotation_summary["compound_nt_var"] = annotation_summary["description"] + annotation_summary["REF"] + annotation_summary["POS"].astype("str") + annotation_summary["ALT"]
    annotation_summary["compound_res_var"] = annotation_summary["description"] + annotation_summary["residues"]
    total_genomic_variants = annotation_summary["compound_nt_var"].nunique()
    total_residue_variants = annotation_summary["compound_res_var"].nunique()

    #for insertions these are annotated as a change on the last ref res !! THINK THIS REGEX NEEDS TO BE USED IN SPEAR ANNOTATION TO CHECK IF REFRES AND ALTRES ARE SAME!!
    annotation_summary['refres'] = annotation_summary["residues"].str.extract('([A-Z])-*[0-9]+-*[a-zA-Z\*\?]+')
    annotation_summary['respos'] = annotation_summary["residues"].str.extract('[A-Z]-*([0-9]+)-*[a-zA-Z\*\?]+')
    annotation_summary['altres'] = annotation_summary["residues"].str.extract('[A-Z]-*[0-9]+-*([a-zA-Z\*\?]+)')

    annotation_summary["respos"] = annotation_summary["respos"].fillna(0).astype(int)
    annotation_summary = annotation_summary.loc[(annotation_summary["refres"] != annotation_summary["altres"]) & (annotation_summary["refres"].isna() == False)]
    baseline_scores = pd.read_csv(f'{args.baseline_scores}', sep = '\t')
    
    orf1ab = ['leader protein', 'nsp2', 'nsp3', 'nsp4', '3C-like proteinase', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 'nsp11', 'RNA-dependent RNA polymerase', 'helicase', "3'-to-5' exonuclease", 'endoRNAse', "2'-O-ribose methyltransferase"]
    structural = ['surface glycoprotein','nucleocapsid phosphoprotein', 'envelope protein', 'membrane glycoprotein']
    accessory = ['ORF3a protein', 'ORF6 protein', 'ORF7a protein', 'ORF7b', 'ORF8 protein', 'ORF10 protein']
    products = orf1ab + structural + accessory
    product_colours = {
        'leader protein': "lightsteelblue", 
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
        'surface glycoprotein': "crimson",
        'nucleocapsid phosphoprotein': "darkred", 
        'envelope protein': "darkred", 
        'membrane glycoprotein': "crimson",
        'ORF3a protein': "cadetblue",
        'ORF6 protein': "darkcyan", 
        'ORF7a protein': "cadetblue", 
        'ORF7b': "darkcyan", 
        'ORF8 protein': "cadetblue",
        'ORF10 protein': "darkcyan"}

    #MAKING A COUNT TABLE OF RAW NUCLEOTIDE VARIANTS - USES THE FILTERED ANNO SUMMARY SO DOES NOT INCLUDE REFRES = ALTRES (BUT SHOULD IT?)
    variants_counts_table = annotation_summary[["sample_id","REF","POS", "ALT","compound_nt_var"]].drop_duplicates(["sample_id", "compound_nt_var"])
    variants_counts_table = variants_counts_table.groupby(["REF","POS", "ALT","compound_nt_var"])[["compound_nt_var"]].count()
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

    buttons = []
    buttons_list = ["Genomic" , "Count"]
    for item in buttons_list:
        if item == "Genomic":
            variants_counts_table = variants_counts_table.sort_values(by = "POS", ascending = True)

        elif item == "Count":
            variants_counts_table = variants_counts_table.sort_values(by = "count", ascending = False)
        buttons.append(dict(
                label = item,
                method = 'restyle',
                args = [
                    {"cells": {
                        "values": [variants_counts_table["REF"], variants_counts_table["POS"], variants_counts_table["ALT"], variants_counts_table["count"]],
                        "fill" : dict(color = 'lavender'),
                        "align":'left'
                    }}]))

    variants_table.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                active = 1,
                bgcolor = "white",
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x = 0.1,
                y = 1.15,
                xanchor="left",
                yanchor="top")
                ])
    variants_table.update_layout(
        annotations=[
            dict(
                text="Sort:", showarrow=False,
                x=0, y=1.1, yref="paper", align="left")
            ]
        )

    variants_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})  
    variants_table_plt = offline.plot(variants_table,output_type='div', include_plotlyjs = True, config = {'displaylogo': False}) #including plotly js with this plot and not with any future ones to keep html size down 

    #MAKING A COUNT TABLE OF RESIDUE CHANGES - USES THE FILTERED ANNO SUMMARY SO DOES NOT INCLUDE REFRES = ALTRES
    residues_counts_table = annotation_summary[["sample_id", "description", "residues", "compound_res_var"]].drop_duplicates(["sample_id", "compound_res_var"])
    residues_counts_table_grouped = residues_counts_table.groupby(["description", "residues"])[["residues"]].count()
    residues_counts_table_grouped.columns = ["count"]
    residues_counts_table_grouped = residues_counts_table_grouped.sort_values("count", axis = 0, ascending = False)
    residues_counts_table_grouped = residues_counts_table_grouped.reset_index()
    residues_counts_table_respos = pd.merge(left = residues_counts_table_grouped, right = annotation_summary[["description", "residues", "respos"]], left_on = ["description", "residues"] , right_on = ["description", "residues"], how = "left")
    residues_counts_table_respos = residues_counts_table_respos.groupby(["description","residues", "respos"]).first().reset_index() #remove duplicates from the merge. 
    residues_counts_table_respos["description"] = residues_counts_table_respos["description"].astype("category")
    residues_counts_table_respos["description"] = residues_counts_table_respos["description"].cat.set_categories(products)
    residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["count", "description", "respos"], ascending = False) #default
    residues_table = go.Figure(data=[go.Table(
        header=dict(values=["Product", "Residue", "Count"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[residues_counts_table_respos["description"], residues_counts_table_respos["residues"], residues_counts_table_respos["count"]],
                   fill_color='lavender',
                   align='left'))
    ])  

    buttons = []
    buttons_list = ["Genomic" , "Count"]
    for item in buttons_list:
        if item == "Genomic":
            residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["description", "respos"])

        elif item == "Count":
            residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["count", "description", "respos"], ascending = False)
        buttons.append(dict(
                label = item,
                method = 'restyle',
                args = [
                    {"cells": {
                        "values": [residues_counts_table_respos["description"], residues_counts_table_respos["residues"], residues_counts_table_respos["count"]],
                        "fill" : dict(color = 'lavender'),
                        "align":'left'
                    }}]))

    residues_table.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                active = 1,
                bgcolor = "white",
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x = 0.05,
                y = 1.15,
                xanchor="left",
                yanchor="top")
                ])

    residues_table.update_layout(
        annotations=[
            dict(
                text="Sort:", showarrow=False,
                x=0, y=1.1, yref="paper", align="left")
            ]
        )
    residues_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
    residues_table_plt = offline.plot(residues_table,output_type='div', include_plotlyjs = False , config = {'displaylogo': False})
    
    #MAKING A SCORES TABLE COLOURED WHERE SAMPLE SCORE SUM EXCEEDS BASELINE 

    scores_cols = baseline_scores.columns.tolist()
    non_displayed_scores = ['total_variants', 'total_residue_variants','consequence_type_variants', 'region_residues', 'domain_residues','ACE2_contact_counts', 'ACE2_contact_score', 'trimer_contact_counts','trimer_contact_score', 'barns_class_variants','bloom_ACE2_max', 'bloom_ACE2_min', 'VDS_max', 'VDS_min', 'serum_escape_max', 'serum_escape_min', 'cm_mAb_escape_all_classes_max','cm_mAb_escape_all_classes_min','mAb_escape_all_classes_max', 'mAb_escape_all_classes_min', 'mAb_escape_class_1_max', 'mAb_escape_class_1_min', 'mAb_escape_class_2_max', 'mAb_escape_class_2_min', 'mAb_escape_class_3_max', 'mAb_escape_class_3_min', 'mAb_escape_class_4_max', 'mAb_escape_class_4_min', 'BEC_RES_max', 'BEC_RES_min', 'BEC_EF_sample']
    displayed_scores_cols = [score for score in scores_cols if score not in non_displayed_scores]
    actual_scores_cols = [score for score in displayed_scores_cols if score != "sample_id"]
    sample_scores = scores_summary[displayed_scores_cols]
    sample_scores = sample_scores.dropna(axis=0, how = "all", subset = actual_scores_cols) #remove empty rows from table to be displayed (do this later for the graph table to allow subtraction of baseline array)
    sample_scores = sample_scores.dropna(axis=1, how = "all") #remove empty cols from table to be displayed (do this later for the graph table to allow subtraction of baseline array)
    sample_scores = pd.concat([baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, displayed_scores_cols], sample_scores])
    sample_scores = sample_scores.reset_index(drop = True)
    #if baseline is also in sample set this can cause problems for the subtraction of a single numpy array from dataframe. Baseline will always be first index value so can fallback on this if multiple matches
    if len(sample_scores.loc[sample_scores["sample_id"] == args.baseline, "sample_id"]) > 1:
        baseline_relative_sample_scores = sample_scores[actual_scores_cols] - sample_scores.loc[0, actual_scores_cols].fillna(0).values.squeeze()
    else:    
        baseline_relative_sample_scores = sample_scores[actual_scores_cols] - sample_scores.loc[sample_scores["sample_id"] == args.baseline, actual_scores_cols].fillna(0).values.squeeze()
    
    sample_scores = sample_scores[sample_scores.index.isin(baseline_relative_sample_scores.index)].reset_index(drop = True)
    baseline_relative_sample_truths = np.greater(baseline_relative_sample_scores.fillna(-1).to_numpy(), 0) #all others should be coloured if greater than baseline
    baseline_relative_sample_colours = np.where(baseline_relative_sample_truths == True, "rgb(250,180,174)", "rgb(179,205,227)")
    baseline_relative_sample_colours[0,:] = np.array(['lavender'] * np.shape(baseline_relative_sample_colours)[1])

    sample_id_colours = np.array([['lavender'] * np.shape(baseline_relative_sample_colours)[0]]).T
    baseline_relative_sample_colours = np.append(sample_id_colours, baseline_relative_sample_colours, axis = 1)
    baseline_relative_sample_colours_df = pd.DataFrame(baseline_relative_sample_colours, columns = displayed_scores_cols)
    sample_scores.set_index("sample_id", drop = False, inplace = True)
    sample_scores.index.name = "index"
    baseline_relative_sample_colours_df.set_index(sample_scores["sample_id"], drop = False, inplace = True)
    sample_scores[actual_scores_cols] = sample_scores[actual_scores_cols].round(2).fillna("").astype("str")

    scores_table = go.Figure(data=[go.Table(
        header=dict(values= [col.replace("_", " ") for col in displayed_scores_cols],
                    fill_color='paleturquoise',
                    align='center'),
        cells=dict(values= [sample_scores[x] for x in displayed_scores_cols],
                fill_color=[baseline_relative_sample_colours_df[x] for x in displayed_scores_cols],
                align='center'))
                ])
        
    buttons = []
    for score in displayed_scores_cols:
        if score == "sample_id":
            asc = True
        else:
            asc = False
        baseline_scores = sample_scores.iloc[[0]]
        baseline_colours = baseline_relative_sample_colours_df.iloc[[0]]
        samples_scores = sample_scores.iloc[1:, :]
        samples_colours = baseline_relative_sample_colours_df.iloc[1:, :]
        samples_scores = samples_scores.replace("", np.nan).sort_values(by = [score], ascending = asc ).replace(np.nan, "")
        samples_colours = samples_colours.reindex(samples_scores.index)
        sorted_scores = pd.concat([baseline_scores, samples_scores])
        sorted_colours = pd.concat([baseline_colours, samples_colours])
        buttons.append(dict(
                label = score,
                method = 'restyle',
                args = [
                    {"cells": {
                        "values": [sorted_scores[x] for x in displayed_scores_cols], 
                        "fill": dict(color = [sorted_colours[x] for x in displayed_scores_cols])}}]))

    scores_table.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                #active = displayed_scores_cols.index("sample_id"),
                direction="down",
                bgcolor = "white",
                pad={"r": 10, "t": 10},
                showactive=True,
                x = 0.025,
                y = 1.15,
                xanchor="left",
                yanchor="top")
                ])
    
    scores_table.update_layout(
        annotations=[
            dict(
                text="Sort:", showarrow=False,
                x=0, y=1.1, yref="paper", align="left")
            ]
        )

    scores_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
    scores_table_plt = offline.plot(scores_table, output_type='div', include_plotlyjs = False , config = {'displaylogo': False})

    if cli_out: #ORDER FOR THIS SHOULD BE DEFAULTED TO CM MAB ESCAPE
        console = Console()
        table = Table(show_header=True, header_style="bold magenta")
        for column in sample_scores.columns:
            table.add_column(column.replace("_"," "), style="dim", width=12)
        cli_baseline_relative_sample_truths = np.where(np.isin(baseline_relative_sample_colours,["rgb(179,205,227)", "lavender"]), False, True)
        for x, y in zip(sample_scores.values, cli_baseline_relative_sample_truths):
            row_value = []
            for (a, b) in zip(x,y):
                colour = "red" if b else "blue"
                value = a if a == a else ""
                row_value.append((f'[bold {colour}]{value}'))
            table.add_row(*row_value)


    #MAKING THE INTERACTIVE PLOTS:
    #use the pkl file of all ORFs to make a heatmap, where you merge the respos from anno df with
    respos_df = pd.read_csv(f'{args.data_dir}/product_mapping.csv')
    orf_boxes = []
    orf_labels = []
    for product in products:
        orf_boxes.append(add_orf_shape(respos_df.loc[respos_df["product"] == product, "position"].min(), respos_df.loc[respos_df["product"] == product, "position"].max(), product, orf1ab, structural, accessory,product_colours))
        orf_labels.append(add_orf_label(respos_df.loc[respos_df["product"] == product, "position"].min(), respos_df.loc[respos_df["product"] == product, "position"].max(), product, orf1ab, structural, accessory))

    count = 0
    if args.product_plots == True:
        table_html_string = []
        table_header = '''<table class="table table-responsive text-center table-hover"><tr><td>Sample ID</td><td>ORF Plot</td><td>Total Variants</td></tr>''' #style="width: auto;"
        table_html_string.append(table_header)
        for sample in track(annotation_summary["sample_id"].unique().tolist(), description="Saving per sample mutation plots..."):
            sample_annotation_summary = annotation_summary.loc[annotation_summary["sample_id"] == sample]
            sample_annotation_summary["respos"] = sample_annotation_summary["respos"].astype(int)
            merged_sample_occurences = pd.merge(sample_annotation_summary, respos_df, left_on = ["description", "respos"], right_on = ["product", "residue-position"], how = "left")
            prods_in_samples = [x for x in products if x in merged_sample_occurences["product"].unique()]
            axis_products = {val:prods_in_samples.index(val) for val in prods_in_samples}
            merged_sample_occurences['axis_pos'] = merged_sample_occurences['product'].map(axis_products)
            merged_sample_occurences = merged_sample_occurences.loc[merged_sample_occurences["residue"].isna() == False]
            merged_sample_occurences["product"] = merged_sample_occurences["product"].astype("category")
            merged_sample_occurences["product"] = merged_sample_occurences["product"].cat.set_categories(products)
            merged_sample_occurences = merged_sample_occurences.sort_values(["product", "respos"])
        
            positive = True
            for prod in prods_in_samples:
                if prod in orf1ab:
                    if positive:
                        merged_sample_occurences.loc[merged_sample_occurences["product"] == prod, "axis_pos"] = np.linspace(start = 8, stop = 1, num = len(merged_sample_occurences.loc[merged_sample_occurences["product"] == prod]))
                        merged_sample_occurences.loc[merged_sample_occurences["product"] == prod, "textpos"] = "top center"
                        positive = False
                    else:
                        merged_sample_occurences.loc[merged_sample_occurences["product"] == prod, "axis_pos"] = np.linspace(start = -8, stop = -1, num = len(merged_sample_occurences.loc[merged_sample_occurences["product"] == prod]))
                        merged_sample_occurences.loc[merged_sample_occurences["product"] == prod, "textpos"] = "bottom center"
                        positive = True

            #rest are direct
            merged_sample_occurences.loc[merged_sample_occurences["product"] == "surface glycoprotein", "axis_pos"] = np.linspace(start = 13, stop = 5, num = len(merged_sample_occurences.loc[merged_sample_occurences["product"] == "surface glycoprotein"]))
            merged_sample_occurences.loc[merged_sample_occurences["product"] == "surface glycoprotein", "textpos"] = "top center"
            merged_sample_occurences.loc[merged_sample_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"]), "axis_pos"] = np.linspace(start = 13, stop = 5, num = len(merged_sample_occurences.loc[merged_sample_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"])]))
            merged_sample_occurences.loc[merged_sample_occurences["product"].isin([x for x in structural if x != "surface glycoprotein"]), "textpos"] = "top center"
            merged_sample_occurences.loc[merged_sample_occurences["product"].isin(accessory), "axis_pos"] = np.linspace(start = -5, stop = -13, num = len(merged_sample_occurences.loc[merged_sample_occurences["product"].isin(accessory)]))
            merged_sample_occurences.loc[merged_sample_occurences["product"].isin(accessory), "textpos"] = "bottom center"
            lines = merged_sample_occurences.apply(lambda x: add_mutation_line(x.position, x.axis_pos, x["product"], orf1ab, structural, accessory),axis = 1).values.tolist()
            combo = orf_boxes + lines
            
            sample_orf_plot = go.Figure()
            sample_orf_plot.add_trace(go.Scatter(
                x=merged_sample_occurences["position"], 
                y=merged_sample_occurences["axis_pos"], 
                text = merged_sample_occurences["residues"],
                mode = "markers+text", 
                textposition= merged_sample_occurences["textpos"]
                ))

            sample_orf_plot.update_traces(
                marker={
                    "size" : 12,
                    "color" : merged_sample_occurences["VDS"],
                    "colorbar" : dict(title = "VDS"),
                    "line" : {"width" : 2 , "color": 'DarkSlateGrey'}
                    },
                selector=dict(mode='markers+text'))
            sample_orf_plot.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',xaxis=dict(showgrid=False, showline = False, zeroline = False),yaxis=dict(showgrid=False, showline = False, zeroline = False)) #title_text='WHAT SHOULD THIS TITLE BE ? ', title_x=0.5, 
            sample_orf_plot.update_layout(
                annotations = orf_labels,
                shapes = combo, 
                title= {"text" : f'{sample} Product Plot'},
                xaxis = {"range" : [0,10000], "fixedrange" : True, "showticklabels" : False, "showgrid" : False},
                yaxis = {"range" : [-14,18], "fixedrange" : True, "showticklabels" : False, "showgrid" : False})

            # Add dropdown
            sample_orf_plot.update_layout(
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
                        bgcolor = "white",
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
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.2,
                        xanchor="left",
                        y=1.08,
                        yanchor="top"
                    ),
                    dict(
                        buttons=list([
                            dict(
                                args=[
                                {"marker" : {
                                        "size" : 12,
                                        "color" : merged_sample_occurences["VDS"],
                                        "colorbar" : dict(title = "VDS"),
                                        "line" : {"width" : 2 , "color": 'DarkSlateGrey'}}}],
                                label="Bloom ACE 2",
                                method="restyle"),
                            dict(
                                args = [{
                                    "marker" : {
                                        "size" : 12,
                                        "color" : merged_sample_occurences["bloom_ace2"],
                                        "colorbar" : dict(title = "Bloom ACE2 Score"),
                                        "line" : {"width" : 2 , "color": 'DarkSlateGrey'}}}],
                                label="Bloom ACE 2",
                                method="restyle"),
                        ]),
                        active = 0,
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.3,
                        xanchor="left",
                        y=1.08,
                        yanchor="top"
                    ),
                    dict(
                        buttons=list([
                            dict(
                                args=[
                                    {"xaxis" : {"range" : [0,10000], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" : [-14,18], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="All",
                                method="relayout"),
                            dict(
                                args = [
                                    {"xaxis" : {"range" :[7000,8400], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[3.8,14], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="Spike Protein",
                                method="relayout"),
                        ]),
                        active = 0,
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.4,
                        xanchor="left",
                        y=1.08,
                        yanchor="top"
                    ),
                ]
            )
            sample_orf_plot.write_html(f'{args.output_dir}/plots/product_plots/{sample}_orf_plot.html', include_plotlyjs = f'../plotly/plotly-2.8.3.min.js')
            sample_total_variants = merged_sample_occurences["residues"].count()
            sample_table_row = f'<tr><td>{sample}</td><td><a href={args.output_dir}/plots/product_plots/{sample}_orf_plot.html>{sample} ORF Plot</a></td><td>TEST</td>'
            table_html_string.append(sample_table_row)
        table_html_string.append("</table>")
        table_html_string = "".join(table_html_string)
        table_card = '''
            <div class="row mt-2">
                <div class="col-12">
                    <div class="card">
                        <div class="card-header" id="headingOne">
                            <h5 class="mb-2">
                                <button class="btn btn-link" data-toggle="collapse" data-target="#collapseOne" aria-expanded="false" aria-controls="collapseOne">
                                ORF Product Plots
                                </button>
                            </h5>
                        </div>
                        <div id="collapseOne" class="collapse show" aria-labelledby="headingOne">
                            <div class="card-body">
                                ''' + table_html_string + '''
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        '''
    else:
        table_card = ""

    #heatmap ________
    
    anno_merge = pd.merge(respos_df, annotation_summary, left_on = ["product", "residue-position"], right_on = ["description", "respos"], how = "left")
    anno_merge.set_index("residues")
    anno_merge.to_csv("~/Desktop/anno_merge.csv")
    scores_cols = ["bloom_ace2", "VDS","serum_escape", "mAb_escape", "cm_mAb_escape", "mAb_escape_class_1", "mAb_escape_class_2", "mAb_escape_class_3", "mAb_escape_class_4", "BEC_RES"]
    scores_z_max = {"bloom_ace2" : 0.3, "VDS": 0.518040609 , "serum_escape" : 1 , "mAb_escape" : 1, "cm_mAb_escape" : 1, "mAb_escape_class_1" : 1, "mAb_escape_class_2" : 1, "mAb_escape_class_3" : 1, "mAb_escape_class_4" : 1, "BEC_RES" : 1}
    scores_z_min = {"bloom_ace2" : -4.84, "VDS" : -0.712636025 ,"serum_escape" : 0 , "mAb_escape" : 0, "cm_mAb_escape" : 0, "mAb_escape_class_1" : 0, "mAb_escape_class_2" : 0, "mAb_escape_class_3" : 0, "mAb_escape_class_4" : 0, "BEC_RES" : 0}
    scores_z_mid = {"serum_escape" : 0.5, "mAb_escape" : 0.5, "cm_mAb_escape" : 0.5, "mAb_escape_class_1" : 0.5, "mAb_escape_class_2" : 0.5, "mAb_escape_class_3" : 0.5, "mAb_escape_class_4" : 0.5, "BEC_RES" : 0.5}
    anno_merge["text_var"] = anno_merge["sample_id"] + ": " + anno_merge["residues"]
    displayed_scores = []
    heatmap = go.Figure()
    heatmap_all = go.Figure()
    for score in scores_cols:
        if anno_merge[score].isna().all() == False:
            displayed_scores.append(score)
            #color bar titles could have a dictionary mapping their underscored titles to a better title, then have an overall plot title ? 
            if score == "VDS":
                heatmap_text = anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"text_var"].values.tolist()
                heatmap_text = [text if text not in [np.nan, "nan"] else "No mutation" for text in heatmap_text]
                heatmap_all_text = anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"text_var"].values.tolist()
                heatmap_all_text = [text if text not in [np.nan, "nan"] else "No mutation" for text in heatmap_all_text]
                heatmap.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"respos"].astype("Int64").astype("str").values.tolist(),
                        'text' : heatmap_text,
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score]
                    }))
                heatmap_all.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913),"residue-position"].astype("Int64").astype("str").values.tolist(),
                        'text' : heatmap_all_text,
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score]
                    }))
            elif score == "bloom_ace2":
                heatmap.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score]
                    }))
                heatmap_all.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "RdBu_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score]
                    }))
            elif score == "cm_mAb_escape":                
                heatmap.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : True,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score],
                        "zmid" : scores_z_mid[score]
                        }))
                heatmap_all.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : True,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score],
                        "zmid" : scores_z_mid[score]
                        }))
            else:
                heatmap.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score],
                        "zmid" : scores_z_mid[score]
                        }))
                heatmap_all.add_trace(go.Heatmap(
                    {
                        'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                        'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"sample_id"].values.tolist(),
                        'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                        'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                        'texttemplate' : "%{text}",
                        'visible' : False,
                        "colorscale" : "hot_r",
                        "colorbar" : dict(title=score), 
                        "zmin" : scores_z_min[score],
                        "zmax" : scores_z_max[score],
                        "zmid" : scores_z_mid[score]
                        }))

    layout = {"title" : dict(text = "Class Masked mAb Escape", x = 0.5), 
        "xaxis" : {"title": "Sample" , "showticklabels" : False, "showgrid" : False},
        "yaxis" : {"title": "Residue Position" ,"tickformat": '.0f', "showgrid" : False}}
    heatmap.update_layout(layout)

    heatmap_all.update_layout(layout)

    trace_list = [True] * len(displayed_scores)
    count = 0
    buttons = []
    for score in displayed_scores:
        trace_list = [False] * len(displayed_scores)
        trace_list[count] = True
        buttons.append(dict(
                        label = score,
                        method = 'update',
                        args = [
                            {'visible': trace_list},
                            {'title': dict(text = score, x = 0.5)},
                            {'yaxes_type' : 'category'}]))
        count += 1
    heatmap.update_layout(paper_bgcolor = 'rgba(0,0,0,0)')  
    heatmap_all.update_layout(paper_bgcolor = 'rgba(0,0,0,0)')
    heatmap.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                active = displayed_scores.index("cm_mAb_escape"),
                bgcolor = "white",
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x = 0.01,
                y = 1.3,
                xanchor="left",
                yanchor="top")
                ])
    heatmap_all.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                active = displayed_scores.index("cm_mAb_escape"),
                bgcolor = "white",
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x = 0.01,
                y = 1.3,
                xanchor="left",
                yanchor="top")
                ])

    heatmap_html = io.to_html(heatmap, full_html = False, include_plotlyjs=False ,config = {'displaylogo': False})
    heatmap_all_html = io.to_html(heatmap, full_html = False, include_plotlyjs=False ,config = {'displaylogo': False})
    heatmap_all.show()
    #maybe a blobbogram of median min max and a baseline line? could have a dropdown on the horizontal bar plot to change the data vis style. 
    
    #################### HTML FORMATTING #####################
    
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <nav class="navbar navbar-expand-md">
	            <div class="container-fluid">	    
		            <div class="collapse navbar-collapse " id="navbarSupportedContent">
		                <ul class="navbar-nav me-auto order-0 ">
			                <li class="nav-item">
			                    <a class="nav-link" href="#"><img src="images/SPEAR_smallest.png" class="img-fluid" alt="Responsive image"></a>
			                </li>
		                </ul>
		                <div class="mx-auto">
			                <a class="navbar-brand text-center" href="#">SPEAR Summary Report</a>
		                </div>
		            </div>
	            </div>
	        </nav>
            <div class="container-fluid">
                <div class = "row mt-2">
                    <div class = "col-12">
                        <div class="card">
                            <div class="card-body">
                            <p class="p1 text-justify">From a total of ''' + str(args.input_samples) + ''' input samples, ''' + str(args.qc_samples) + ''' passed QC and were annotated by SPEAR. A total of ''' + str(total_genomic_variants) + ''' nucleotide variants were identified across all samples, which resulted in  ''' + str(total_residue_variants) + ''' AA missense, deletion or insertions which are listed and evaluated below. </p>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="row mt-2">
                    <div class="col-6">
                        <div class="card">
                            <div class="card-header">Nucleotide Variants</div>
                            <div class="card-body">
                                <div>
                                ''' + variants_table_plt + '''
                                </div>
                            </div>
                            <div class="card-footer">TEST TEXT</div>
                        </div>
                    </div>
                    <div class="col-6">
                        <div class="card">
                            <div class="card-header">Amino Acid Variants</div>
                            <div class="card-body">
                                <div>
                                ''' + residues_table_plt + '''
                                </div>
                            </div>
                            <div class="card-footer">TEST TEXT</div>
                        </div>
                    </div>
                </div>
                <div class = "row mt-2">
                    <div class="col-12">
                        <div class="card">
                            <div class="card-header">Residue Score Heatmap</div>
                            <div class="card-body">
                                ''' + heatmap_html + '''
                            </div>
                            <div class="card-footer">
                                The above plot shows a heatmap of all residue positions for which the residue was changed in at least one sample, and the associated score. <a href="https://github.com/m-crown/SPEAR">Explanation of scores</a> <a href="https://github.com/m-crown/SPEAR">Full screen</a>
                            </div>
                        </div>
                    </div>
                </div>
                <div class = "row mt-2">
                     <div class="col-12">
                        <div class="card">
                            <div class="card-header">Per Sample Scores Summary</div>
                            <div class="card-body">
                                ''' + scores_table_plt + '''
                            </div>
                            </div>
                        </div>
                    </div> 
                </div>       
                ''' + table_card + '''
            </div>
            <!-- Footer -->
            <footer class="page-footer font-small teal pt-4">
                <!-- Footer Text -->
                <div class="container-fluid text-center text-md-left">
                    <!-- Grid row -->
                    <div class="row">
                    <!-- Grid column -->
                    <div class="col-md-6 mt-md-0 mt-3">
                        <!-- Content -->
                        <h5 class="text-uppercase font-weight-bold">Footer text 1</h5>
                        <p>This report was generated using SPEAR using v0.4.</p>
                    </div>
                    <!-- Grid column -->
                    <hr class="clearfix w-100 d-md-none pb-3">

                    <!-- Grid column -->
                    <div class="col-md-6 mb-md-0 mb-3">

                        <!-- Content -->
                        <h5 class="text-uppercase font-weight-bold">Footer text 2</h5>
                        <p>Lorem ipsum, dolor sit amet consectetur adipisicing elit. Optio deserunt fuga perferendis modi earum
                        commodi aperiam temporibus quod nulla nesciunt aliquid debitis ullam omnis quos ipsam, aspernatur id
                        excepturi hic.</p>

                    </div>
                    <!-- Grid column -->

                    </div>
                    <!-- Grid row -->

                    </div>
                    <!-- Footer Text -->

                </footer>
                <!-- Footer -->
 
        <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
        <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
        </body>
    </html>'''
    
    copyfile(f'{args.images_dir}/SPEAR_smallest.png', f'{args.output_dir}/images/SPEAR_smallest.png')
    copyfile(f'{args.images_dir}/GitHub-Mark-Light-64px.png', f'{args.output_dir}/images/GitHub-Mark-Light-32px.png')
    f = open(f'{args.output_dir}/report.html','w')
    f.write(html_string)
    f.close()

    #after all processing is complete, print the table of results. 
    #console.print(table)
if __name__ == "__main__":
    main()