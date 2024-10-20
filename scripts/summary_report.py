#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objects as go
import plotly.offline as offline
import argparse
from pathlib import Path
import numpy as np
from shutil import copyfile
from numpy.random import seed
from rich.console import Console
from rich.table import Table
from rich.progress import track
import os
import csv

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


def summary_report(args):

    seed(42069)
    Path(f'{args.output_dir}/images').mkdir(parents=True, exist_ok=True)
    Path(f'{args.output_dir}/plots/plotly').mkdir(parents=True, exist_ok=True)
    Path(f'{args.output_dir}/plots/product_plots').mkdir(parents=True, exist_ok=True)
    copyfile(f'{args.scripts_dir}/plotly-2.8.3.min.js', f'{args.output_dir}/plots/plotly/plotly-2.8.3.min.js')
    
    report_date = pd.to_datetime('today').strftime('%Y-%m-%d')

    console = Console()

    scores_summary = pd.read_csv(f'{args.score_summary}', sep = '\t')
    annotation_summary = pd.read_csv(f'{args.annotation_summary}', sep = '\t')
    try:
        lineages = pd.read_csv(f'{args.pangolin_report}', sep = ",")
        run_pangolin = True
        if len(lineages) == 0:
            lineages = pd.DataFrame(data = {"taxon" : scores_summary["sample_id"], "lineage" : ""}) 
    except pd.errors.EmptyDataError:
        lineages = pd.DataFrame(data = {"taxon" : scores_summary["sample_id"], "lineage" : ""})
        run_pangolin = False

    with open(f'{args.pangolin_command}') as file:
        pangolin_command = [line.rstrip() for line in file]
    pangolin_params = pangolin_command[0]
    pangolin_versions = ' '.join(pangolin_command[1:])

    with open(args.spear_params) as file:
        spear_params = file.readline().rstrip()

    spear_params_list = spear_params.split(',')
    spear_params_dict = {param.split(':')[0]: param.split(':')[1] for param in spear_params_list}
    global_n = float(spear_params_dict["global_n"])
    s_n = float(spear_params_dict["s_n"])
    s_contig = int(spear_params_dict["s_contig"])
    rbd_n = int(spear_params_dict["rbd_n"])
    
    spear_qc_info = pd.read_csv(f"{args.spear_qc_info}", sep = "\t")
    spear_version = spear_qc_info.spear_version.values[0]
    qc_samples = spear_qc_info.passing_samples.astype(int).values[0]

    annotation_summary["compound_nt_var"] = annotation_summary["description"] + annotation_summary["REF"] + annotation_summary["POS"].astype("str") + annotation_summary["ALT"]
    annotation_summary["compound_res_var"] = annotation_summary["description"] + annotation_summary["residues"]
    total_genomic_variants = annotation_summary["compound_nt_var"].nunique()

    #for insertions these are annotated as a change on the last ref res !! THINK THIS REGEX NEEDS TO BE USED IN SPEAR ANNOTATION TO CHECK IF REFRES AND ALTRES ARE SAME!!
    annotation_summary['refres'] = annotation_summary["residues"].str.extract('([A-Z\*])-*[0-9]+-*[a-zA-Z\*\?]+')
    annotation_summary['respos'] = annotation_summary["residues"].str.extract('[A-Z\*]-*([0-9]+)-*[a-zA-Z\*\?]+')
    annotation_summary['altres'] = annotation_summary["residues"].str.extract('[A-Z\*]-*[0-9]+-*([a-zA-Z\*\?]+)')

    #MAKING A COUNT TABLE OF RAW NUCLEOTIDE VARIANTS
    variants_counts_table = annotation_summary[["sample_id","REF","POS", "ALT","compound_nt_var"]].drop_duplicates(["sample_id", "compound_nt_var"])
    variants_counts_table = variants_counts_table.groupby(["REF","POS", "ALT","compound_nt_var"])[["compound_nt_var"]].count()
    variants_counts_table.columns = ["count"]
    variants_counts_table = variants_counts_table.sort_values("count", axis = 0 , ascending = False)
    variants_counts_table = variants_counts_table.reset_index()
    variants_counts_table["percentage_samples"] = (variants_counts_table["count"] / int(qc_samples)) * 100
    variants_counts_table["percentage_samples"] = variants_counts_table["percentage_samples"].round(2)
     
    variants_table = go.Figure(data=[go.Table(
        header=dict(values= ["POS","REF", "ALT", "count", "% Samples"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[variants_counts_table["REF"], variants_counts_table["POS"], variants_counts_table["ALT"], variants_counts_table["count"], variants_counts_table["percentage_samples"]],
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
                        "values": [variants_counts_table["REF"], variants_counts_table["POS"], variants_counts_table["ALT"], variants_counts_table["count"], variants_counts_table["percentage_samples"]],
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
    #variants_table_plt = offline.plot(variants_table,output_type='div', include_plotlyjs = f'plots/plotly/plotly-2.8.3.min.js', config = {'displaylogo': False}) #including plotly js with this plot and not with any future ones to keep html size down 
    variants_table.write_html(f'{args.output_dir}/plots/nt_variants_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')

    annotation_summary["respos"] = annotation_summary["respos"].fillna(0).astype(int)
    annotation_summary = annotation_summary.loc[(annotation_summary["refres"] != annotation_summary["altres"]) & (annotation_summary["refres"].isna() == False)]
    total_residue_variants = annotation_summary["compound_res_var"].nunique()
    
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

    product_order = [
        'leader protein', 
        'nsp2', 
        'nsp3', 
        'nsp4', 
        '3C-like proteinase', 
        'nsp6', 
        'nsp7', 
        'nsp8', 
        'nsp9', 
        'nsp10', 
        'nsp11', 
        'RNA-dependent RNA polymerase', 
        'helicase', 
        "3'-to-5' exonuclease", 
        'endoRNAse', 
        "2'-O-ribose methyltransferase",
        'surface glycoprotein',
        'ORF3a protein',
        'envelope protein', 
        'membrane glycoprotein',
        'ORF6 protein', 
        'ORF7a protein', 
        'ORF7b', 
        'ORF8 protein',
        'ORF10 protein',
        'nucleocapsid phosphoprotein']

    #MAKING A COUNT TABLE OF RESIDUE CHANGES - USES THE FILTERED ANNO SUMMARY SO DOES NOT INCLUDE REFRES = ALTRES
    residues_counts_table = annotation_summary[["sample_id", "description", "residues", "compound_nt_var", "compound_res_var"]].drop_duplicates(["sample_id", "compound_res_var"])
    residues_counts_table_grouped = residues_counts_table.groupby(["compound_nt_var", "description", "residues"])[["residues"]].count()
    residues_counts_table_grouped.columns = ["count"]
    residues_counts_table_grouped = residues_counts_table_grouped.sort_values("count", axis = 0, ascending = False)
    residues_counts_table_grouped = residues_counts_table_grouped.reset_index()
    residues_counts_table_respos = pd.merge(left = residues_counts_table_grouped, right = annotation_summary[["POS", "REF", "ALT", "compound_nt_var", "description", "residues", "respos"]], left_on = ["compound_nt_var", "description", "residues"] , right_on = ["compound_nt_var", "description", "residues"], how = "left")
    residues_counts_table_respos = residues_counts_table_respos.groupby(["compound_nt_var", "description","residues", "respos"]).first().reset_index() #remove duplicates from the merge. 
    residues_counts_table_respos["description"] = residues_counts_table_respos["description"].astype("category")
    residues_counts_table_respos["description"] = residues_counts_table_respos["description"].cat.set_categories(products)
    residues_counts_table_respos["order"] = residues_counts_table_respos["description"].apply(lambda x: product_order.index(x))
    residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["count", "order", "respos"], ascending = [False, True, True]) #default
    residues_counts_table_respos["percentage_samples"] = (residues_counts_table_respos["count"] / int(qc_samples)) * 100
    residues_counts_table_respos["percentage_samples"] = residues_counts_table_respos["percentage_samples"].round(2)

    residues_table = go.Figure(data=[go.Table(
        header=dict(values=["POS", "REF", "ALT", "Product", "Residue", "Count", "% Samples"],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[residues_counts_table_respos["POS"],residues_counts_table_respos["REF"],residues_counts_table_respos["ALT"], residues_counts_table_respos["description"], residues_counts_table_respos["residues"], residues_counts_table_respos["count"], residues_counts_table_respos["percentage_samples"]],
                   fill_color='lavender',
                   align='left'))
    ])  

    buttons = []
    buttons_list = ["Genomic" , "Count"]
    for item in buttons_list:
        if item == "Genomic":
            residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["order", "respos"])

        elif item == "Count":
            residues_counts_table_respos = residues_counts_table_respos.sort_values(by = ["count", "order", "respos"], ascending = [False, True, True])
        buttons.append(dict(
                label = item,
                method = 'restyle',
                args = [
                    {"cells": {
                        "values": [residues_counts_table_respos["POS"],residues_counts_table_respos["REF"],residues_counts_table_respos["ALT"],residues_counts_table_respos["description"], residues_counts_table_respos["residues"], residues_counts_table_respos["count"], residues_counts_table_respos["percentage_samples"]],
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
    residues_table_plt = offline.plot(residues_table,output_type='div', include_plotlyjs = f'plots/plotly/plotly-2.8.3.min.js', config = {'displaylogo': False})

    #making a domain counts table
    if annotation_summary["domain"].replace("", np.nan).isnull().all() == False: 

        domain_counts_table_all = annotation_summary.loc[annotation_summary["domain"].replace("", np.nan).isnull() == False , ["description", "residues", "respos", "domain"]].copy()
        domain_counts_table_all["domain"] = domain_counts_table_all["domain"].replace("_", " ", regex = True)
        domain_counts_table_all_grouped = domain_counts_table_all.groupby(["description", "residues", "respos", "domain"])[["domain"]].count()
        domain_counts_table_all_grouped.columns = ["count"]
        domain_counts_table_all_grouped = domain_counts_table_all_grouped.sort_values("count", axis = 0, ascending = False)
        domain_counts_table_all_grouped = domain_counts_table_all_grouped.reset_index()
        domain_counts_table_all_grouped.rename(columns = {"count" : "total_mutations"}, inplace = True)
        final_domain_counts = domain_counts_table_all_grouped
        final_domain_counts["order"] = final_domain_counts["description"].apply(lambda x: product_order.index(x))
        final_domain_counts = final_domain_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
        domain_table = go.Figure(data=[go.Table(
            header=dict(values=["Product", "Mutation", "Domain", "Samples"],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[final_domain_counts["description"], final_domain_counts["residues"], final_domain_counts["domain"], final_domain_counts["total_mutations"]],
                        fill_color='lavender',
                        align='left'))
        ])  

        buttons = []
        buttons_list = ["Total" , "Product"]
        for item in buttons_list:
            if item == "Total":
                final_domain_counts = final_domain_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
            elif item == "Product": 
                final_domain_counts = final_domain_counts.sort_values(by = ["order", "respos"], ascending = [True, True])
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [final_domain_counts["description"], final_domain_counts["residues"], final_domain_counts["domain"], final_domain_counts["total_mutations"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        domain_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        domain_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        domain_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        domain_table_plt = offline.plot(domain_table,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})
        all_samples_domain = annotation_summary.loc[annotation_summary["domain"].replace("", np.nan).isnull() == False , ["sample_id", "description", "residues", "respos", "domain"]].copy()
        all_samples_domain["domain"] = all_samples_domain["domain"].replace("_", " ", regex = True)
        all_samples_domain["order"] = all_samples_domain["description"].apply(lambda x: product_order.index(x))
        all_samples_domain.sort_values(by = ["sample_id", "order", "respos"], ascending = [True, True, True],inplace = True)
        all_samples_domain_table = go.Figure(data=[go.Table(
                header=dict(values=["Sample", "Product", "Domain", "residue"],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[all_samples_domain["sample_id"],all_samples_domain["description"],all_samples_domain["domain"], all_samples_domain["residues"]],
                            fill_color='lavender',
                            align='left'))
            ])
        buttons = []
        buttons_list = ["Sample" , "Product"]
        for item in buttons_list:
            if item == "Sample":
                all_samples_domain = all_samples_domain.sort_values(by = ["sample_id", "order", "respos"], ascending = [True, True, True])
            elif item == "Product": 
                all_samples_domain = all_samples_domain.sort_values(by = ["order", "respos"], ascending = [True, True])
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [all_samples_domain["sample_id"],all_samples_domain["description"],all_samples_domain["domain"], all_samples_domain["residues"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        all_samples_domain_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        all_samples_domain_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        all_samples_domain_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        all_samples_domain_table.write_html(f'{args.output_dir}/plots/samples_domain_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        domain_message = '''Table of mutation counts in protein domains. Total mutations: the total number of mutations in this domain seen across all samples. Unique mutations: total number of unique residue mutations observed across all samples.
                            A full table of protein domain mutation in each sample can be found <a href="plots/samples_domain_table.html">here</a> and in <code>spear_annotation_summary.tsv</code>'''

    else:
        domain_table_plt = "No domain mutations present in any samples."
        domain_message = ""

    #making a feature counts table
    if annotation_summary["feature"].replace("", np.nan).isnull().all() == False: 

        feature_counts_table_all = annotation_summary.loc[annotation_summary["feature"].replace("", np.nan).isnull() == False , ["description", "residues", "respos", "domain", "feature"]].copy()
        feature_counts_table_all[["feature", "domain"]] = feature_counts_table_all[["feature", "domain"]].replace({" " : ",", np.nan : ""}, regex = True)
        feature_counts_table_all[["feature", "domain"]] = feature_counts_table_all[["feature", "domain"]].replace({"_": " "}, regex = True)
        feature_counts_table_all_grouped = feature_counts_table_all.groupby(["description", "domain", "residues", "respos", "feature"])[["feature"]].count()
        feature_counts_table_all_grouped.columns = ["count"]
        feature_counts_table_all_grouped = feature_counts_table_all_grouped.sort_values("count", axis = 0, ascending = False)
        feature_counts_table_all_grouped = feature_counts_table_all_grouped.reset_index()
        feature_counts_table_all_grouped.rename(columns = {"count" : "total_mutations"}, inplace = True)
        final_feature_counts = feature_counts_table_all_grouped
        final_feature_counts["order"] = final_feature_counts["description"].apply(lambda x: product_order.index(x))
        final_feature_counts = final_feature_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
        feature_table = go.Figure(data=[go.Table(
            header=dict(values=["Product", "Mutation", "Domain", "Feature", "Samples"],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[final_feature_counts["description"], final_feature_counts["residues"], final_feature_counts["domain"],final_feature_counts["feature"], final_feature_counts["total_mutations"]],
                        fill_color='lavender',
                        align='left'))
        ])  

        buttons = []
        buttons_list = ["Total", "Product"]
        for item in buttons_list:
            if item == "Total":
                final_feature_counts = final_feature_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
            elif item == "Product": 
                final_feature_counts = final_feature_counts.sort_values(by = ["order", "respos"], ascending = [True, True])
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [final_feature_counts["description"], final_feature_counts["residues"], final_feature_counts["domain"],final_feature_counts["feature"], final_feature_counts["total_mutations"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        feature_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        feature_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        feature_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        feature_table_plt = offline.plot(feature_table,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})
        all_samples_feature = annotation_summary.loc[annotation_summary["feature"].replace("", np.nan).isnull() == False , ["sample_id", "description", "residues", "respos", "domain", "feature"]].copy()
        all_samples_feature[["feature", "domain"]] = all_samples_feature[["feature", "domain"]].replace({" " : ",", np.nan : ""}, regex = True)
        all_samples_feature[["feature", "domain"]] = all_samples_feature[["feature", "domain"]].replace({"_": " "}, regex = True)
        all_samples_feature["order"] = all_samples_feature["description"].apply(lambda x: product_order.index(x))
        all_samples_feature.sort_values(by = ["sample_id", "order", "respos"], inplace = True)
        all_samples_feature_table = go.Figure(data=[go.Table(
                header=dict(values=["Sample", "Product", "Domain", "Feature", "Residue"],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[all_samples_feature["sample_id"],all_samples_feature["description"],all_samples_feature["domain"],all_samples_feature["feature"], all_samples_feature["residues"]],
                            fill_color='lavender',
                            align='left'))
            ])
        buttons = []
        buttons_list = ["Sample" , "Product"]
        for item in buttons_list:
            if item == "Sample":
                all_samples_feature = all_samples_feature.sort_values(by = ["sample_id", "order", "respos"], ascending = True)
            elif item == "Product": 
                all_samples_feature = all_samples_feature.sort_values(by = ["order", "respos"], ascending = True)
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [all_samples_feature["sample_id"], all_samples_feature["description"],all_samples_feature["domain"], all_samples_feature["feature"],all_samples_feature["residues"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        all_samples_feature_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        all_samples_feature_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        all_samples_feature_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        all_samples_feature_table.write_html(f'{args.output_dir}/plots/samples_feature_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        feature_message = '''Table of mutation counts of protein features (in domains). Total mutations: the total number of mutations in this feature seen across all samples. Unique mutations: total number of unique residue mutations observed across all samples in this feature.
                            A full table of protein feature mutation in each sample can be found <a href="plots/samples_feature_table.html">here</a> and in <code>spear_annotation_summary.tsv</code>'''

    else:
        final_feature_counts = pd.DataFrame({"sample_id" : [] , "description" : [], "respos" : []})
        feature_table_plt = "No domain mutations present in any samples."
        feature_message = ""
    
    #making a furin cleavage site specific table
    if len(final_feature_counts.loc[(final_feature_counts["description"] == "surface glycoprotein") & (final_feature_counts["respos"] >= 675) & (final_feature_counts["respos"] <= 695)]) != 0:

        furin_feature_counts = final_feature_counts.loc[(final_feature_counts["description"] == "surface glycoprotein") & (final_feature_counts["respos"] >= 675) & (final_feature_counts["respos"] <= 695)].copy()
        furin_table = go.Figure(data=[go.Table(
            header=dict(values=["Product", "Mutation", "Domain", "Feature", "Samples"],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[furin_feature_counts["description"], furin_feature_counts["residues"], furin_feature_counts["domain"],furin_feature_counts["feature"], furin_feature_counts["total_mutations"]],
                        fill_color='lavender',
                        align='left'))
        ])  

        buttons = []
        buttons_list = ["Total", "Product"]
        for item in buttons_list:
            if item == "Total":
                furin_feature_counts = furin_feature_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
            elif item == "Product": 
                furin_feature_counts = furin_feature_counts.sort_values(by = ["order", "respos"], ascending = [True, True])
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [furin_feature_counts["description"], furin_feature_counts["residues"], furin_feature_counts["domain"],furin_feature_counts["feature"], furin_feature_counts["total_mutations"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        furin_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        furin_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        furin_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        furin_table_plt = offline.plot(furin_table,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})
        all_samples_furin = all_samples_feature[(all_samples_feature["description"] == "surface glycoprotein") & (all_samples_feature["respos"] >= 680) & (all_samples_feature["respos"] <= 690)].copy()
        all_samples_furin[["feature", "domain"]] = all_samples_furin[["feature", "domain"]].replace({" " : ",", np.nan : ""}, regex = True)
        all_samples_furin[["feature", "domain"]] = all_samples_furin[["feature", "domain"]].replace({"_": " "}, regex = True)
        all_samples_furin["order"] = all_samples_furin["description"].apply(lambda x: product_order.index(x))
        all_samples_furin.sort_values(by = ["sample_id", "order", "respos"], inplace = True)
        all_samples_furin_table = go.Figure(data=[go.Table(
                header=dict(values=["Sample", "Product", "Domain", "Feature", "Residue"],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[all_samples_furin["sample_id"],all_samples_furin["description"],all_samples_furin["domain"],all_samples_furin["feature"],all_samples_furin["residues"]],
                            fill_color='lavender',
                            align='left'))
            ])
        buttons = []
        buttons_list = ["Sample" , "Product"]
        for item in buttons_list:
            if item == "Sample":
                all_samples_furin = all_samples_furin.sort_values(by = ["sample_id", "order", "respos"], ascending = True)
            elif item == "Product": 
                all_samples_furin = all_samples_furin.sort_values(by = ["order", "respos"], ascending = True)
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [all_samples_furin["sample_id"],all_samples_furin["description"],all_samples_furin["domain"],all_samples_furin["feature"],all_samples_furin["residues"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        all_samples_furin_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        all_samples_furin_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        all_samples_furin_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        all_samples_furin_table.write_html(f'{args.output_dir}/plots/samples_furin_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        furin_message = '''Subset of feature table containing mutations in Spike residues 675-695, the region of the S1/S2 furin cleavage site.
                            A full table of furin site mutations in each sample can be found <a href="plots/samples_furin_table.html">here</a> and in <code>spear_annotation_summary.tsv</code>'''

    else:
        furin_table_plt = "No furin mutations present in any samples."
        furin_message = ""

    #making an NTD supersite loops specific table
    if len(annotation_summary.loc[annotation_summary["domain"] == "NTD"]) != 0:
        ntd_counts_table_all = annotation_summary.loc[annotation_summary["domain"] == "NTD" , ["description", "residues", "respos", "domain", "feature"]].copy()
        ntd_counts_table_all[["feature", "domain"]] = ntd_counts_table_all[["feature", "domain"]].replace({" " : ",", np.nan : ""}, regex = True)
        ntd_counts_table_all[["feature", "domain"]] = ntd_counts_table_all[["feature", "domain"]].replace({"_" : " "}, regex = True)
        ntd_counts_table_all_grouped = ntd_counts_table_all.groupby(["description", "domain", "residues", "respos", "feature"])[["feature"]].count()
        ntd_counts_table_all_grouped.columns = ["count"]
        ntd_counts_table_all_grouped = ntd_counts_table_all_grouped.sort_values("count", axis = 0, ascending = False)
        ntd_counts_table_all_grouped = ntd_counts_table_all_grouped.reset_index()
        ntd_counts_table_all_grouped.rename(columns = {"count" : "total_mutations"}, inplace = True)
        final_ntd_counts = ntd_counts_table_all_grouped
        final_ntd_counts["order"] = final_ntd_counts["description"].apply(lambda x: product_order.index(x))
        final_ntd_counts = final_ntd_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
        ntd_table = go.Figure(data=[go.Table(
            header=dict(values=["Product", "Mutation", "Domain", "Feature", "Samples"],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[final_ntd_counts["description"], final_ntd_counts["residues"], final_ntd_counts["domain"],final_ntd_counts["feature"], final_ntd_counts["total_mutations"]],
                        fill_color='lavender',
                        align='left'))
        ])  

        buttons = []
        buttons_list = ["Total", "Product"]
        for item in buttons_list:
            if item == "Total":
                final_ntd_counts = final_ntd_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
            elif item == "Product": 
                final_ntd_counts = final_ntd_counts.sort_values(by = ["order", "respos"], ascending = [True, True])
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [final_ntd_counts["description"], final_ntd_counts["residues"], final_ntd_counts["domain"],final_ntd_counts["feature"], final_ntd_counts["total_mutations"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        ntd_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        ntd_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        ntd_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        ntd_table_plt = offline.plot(ntd_table,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})
        
        all_samples_ntd = annotation_summary.loc[annotation_summary["domain"] == "NTD" , ["sample_id", "description", "residues", "respos", "domain", "feature"]].copy()
        all_samples_ntd[["feature", "domain"]] = all_samples_ntd[["feature", "domain"]].replace({" " : ",", np.nan : ""}, regex = True)
        all_samples_ntd[["feature", "domain"]] = all_samples_ntd[["feature", "domain"]].replace({"_": " "}, regex = True)
        all_samples_ntd["order"] = all_samples_ntd["description"].apply(lambda x: product_order.index(x))
        all_samples_ntd.sort_values(by = ["sample_id", "order", "respos"], inplace = True)
        all_samples_ntd_table = go.Figure(data=[go.Table(
                header=dict(values=["Sample", "Product", "Domain", "Feature", "Residue"],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[all_samples_ntd["sample_id"],all_samples_ntd["description"],all_samples_ntd["domain"],all_samples_ntd["feature"],all_samples_ntd["residues"]],
                            fill_color='lavender',
                            align='left'))
            ])
        buttons = []
        buttons_list = ["Sample" , "Product"]
        for item in buttons_list:
            if item == "Sample":
                all_samples_ntd = all_samples_ntd.sort_values(by = ["sample_id", "order", "respos"], ascending = True)
            elif item == "Product": 
                all_samples_ntd = all_samples_ntd.sort_values(by = ["order", "respos"], ascending = True)
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [all_samples_ntd["sample_id"],all_samples_ntd["description"],all_samples_ntd["domain"],all_samples_ntd["feature"],all_samples_ntd["residues"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        all_samples_ntd_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        all_samples_ntd_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        all_samples_ntd_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        all_samples_ntd_table.write_html(f'{args.output_dir}/plots/samples_ntd_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        ntd_message = '''Mutations in Spike NTD.
                            A full table of NTD mutations in each sample can be found <a href="plots/samples_ntd_table.html">here</a> and in <code>spear_annotation_summary.tsv</code>'''

    else:
        ntd_table_plt = "No NTD mutations present in any samples."
        ntd_message = ""
    
    #making a contacts counts table
    if annotation_summary["contact_type"].replace("", np.nan).isnull().all() == False: 
        contact_counts_table_all = annotation_summary.loc[annotation_summary["contact_type"].replace("", np.nan).isnull() == False , ["description", "residues", "respos", "contact_type"]].copy()
        contact_counts_table_all["contact_type"] = contact_counts_table_all["contact_type"].str.split(" ")
        contact_counts_table_all = contact_counts_table_all.explode("contact_type")
        contact_counts_table_all[["contact", "contact_type"]] =  contact_counts_table_all["contact_type"].str.split(":", expand = True)
        contact_counts_table_all["contact_type"] = contact_counts_table_all["contact_type"].str.split("+")
        contact_counts_table_all = contact_counts_table_all.explode("contact_type")
        contact_counts_table_all["contact"] = contact_counts_table_all["contact"].replace("_", " ", regex = True)
        contact_counts_table_all.loc[contact_counts_table_all["contact_type"].fillna("").str.match(r"[A-Za-z]+$"), "contact_type"] = contact_counts_table_all.loc[contact_counts_table_all["contact_type"].fillna("").str.match(r"[A-Za-z]+$"), "contact_type"] + "_" 
        contact_counts_table_all[["contact_type", "contacting_residues"]] = contact_counts_table_all["contact_type"].str.split("_", n=1, expand = True)
        contact_counts_table_all["contacting_residues"].fillna("", inplace = True)
        contact_counts_table_all["contacting_residues"] = contact_counts_table_all["contacting_residues"].str.replace("_", ",")
        contact_counts_table_all_grouped = contact_counts_table_all.groupby(["description", "residues", "respos", "contact", "contact_type", "contacting_residues"])[["contact_type"]].count()
        contact_counts_table_all_grouped.columns = ["count"]
        contact_counts_table_all_grouped = contact_counts_table_all_grouped.sort_values("count", axis = 0, ascending = False)
        contact_counts_table_all_grouped = contact_counts_table_all_grouped.reset_index()

        contact_counts_table_all_grouped.rename(columns = {"count" : "total_mutations"}, inplace = True)
        final_contact_counts = contact_counts_table_all_grouped
        final_contact_counts["order"] = final_contact_counts["description"].apply(lambda x: product_order.index(x))
        final_contact_counts = final_contact_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
        contacts_table = go.Figure(data=[go.Table(
            header=dict(values=["Product", "Mutation", "Contact", "Contact Type","Contacting Residue(s)", "Samples"],
                        fill_color='paleturquoise',
                        align='left'),
            cells=dict(values=[final_contact_counts["description"],final_contact_counts["residues"],final_contact_counts["contact"],final_contact_counts["contact_type"],final_contact_counts["contacting_residues"], final_contact_counts["total_mutations"]],
                        fill_color='lavender',
                        align='left'))
        ])  

        buttons = []
        buttons_list = ["Total" , "Product"]
        for item in buttons_list:
            if item == "Total":
                final_contact_counts = final_contact_counts.sort_values(by = ["total_mutations", "order", "respos"], ascending = [False, True, True])
            elif item == "Product": 
                final_contact_counts = final_contact_counts.sort_values(by = ["order", "respos"], ascending = True)
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [final_contact_counts["description"],final_contact_counts["residues"],final_contact_counts["contact"],final_contact_counts["contact_type"],final_contact_counts["contacting_residues"], final_contact_counts["total_mutations"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        contacts_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        contacts_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        contacts_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        contacts_table_plt = offline.plot(contacts_table,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})

        all_samples_contacts = annotation_summary.loc[annotation_summary["contact_type"].replace("", np.nan).isnull() == False , ["sample_id", "description", "residues", "respos", "contact_type"]].copy()
        all_samples_contacts["contact_type"] = all_samples_contacts["contact_type"].str.split(" ")
        all_samples_contacts = all_samples_contacts.explode("contact_type")
        all_samples_contacts[["contact", "contact_type"]] =  all_samples_contacts["contact_type"].str.split(":", expand = True)
        all_samples_contacts["contact_type"] = all_samples_contacts["contact_type"].str.split("+")
        all_samples_contacts = all_samples_contacts.explode("contact_type")
        all_samples_contacts["contact"] = all_samples_contacts["contact"].replace("_", " ", regex = True)
        all_samples_contacts.loc[all_samples_contacts["contact_type"].fillna("").str.match(r"[A-Za-z]+$"), "contact_type"] = all_samples_contacts.loc[all_samples_contacts["contact_type"].fillna("").str.match(r"[A-Za-z]+$"), "contact_type"] + "_" 
        all_samples_contacts[["contact_type", "contacting_residues"]] = all_samples_contacts["contact_type"].str.split("_", n=1, expand = True)
        all_samples_contacts["contacting_residues"].fillna("", inplace = True)
        all_samples_contacts["contacting_residues"] = all_samples_contacts["contacting_residues"].str.replace("_", ",")
        all_samples_contacts["order"] = all_samples_contacts["description"].apply(lambda x: product_order.index(x))
        all_samples_contacts.sort_values(by = ["sample_id", "order", "respos"], inplace = True)
        all_samples_contacts_table = go.Figure(data=[go.Table(
                header=dict(values=["Sample", "Product", "Mutation", "Contact", "Contact Type","Contacting Residue(s)"],
                            fill_color='paleturquoise',
                            align='left'),
                cells=dict(values=[all_samples_contacts["sample_id"],all_samples_contacts["description"],all_samples_contacts["residues"],all_samples_contacts["contact"],all_samples_contacts["contact_type"], all_samples_contacts["contacting_residues"]],
                            fill_color='lavender',
                            align='left'))
            ])
        buttons = []
        buttons_list = ["Sample" , "Product"]
        for item in buttons_list:
            if item == "Sample":
                all_samples_contacts = all_samples_contacts.sort_values(by = ["sample_id", "order"], ascending = True)
            elif item == "Product": 
                all_samples_contacts = all_samples_contacts.sort_values(by = "order", ascending = True)
            buttons.append(dict(
                    label = item,
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [all_samples_contacts["sample_id"],all_samples_contacts["description"],all_samples_contacts["residues"],all_samples_contacts["contact"],all_samples_contacts["contact_type"], all_samples_contacts["contacting_residues"]],
                            "fill" : dict(color = 'lavender'),
                            "align":'left'
                        }}]))

        all_samples_contacts_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = 0,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.05,
                    y = 1.15,
                    xanchor="left",
                    yanchor="top")
                    ])

        all_samples_contacts_table.update_layout(
            annotations=[
                dict(
                    text="Sort:", showarrow=False,
                    x=0, y=1.1, yref="paper", align="left")
                ]
            )
        all_samples_contacts_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5)})
        all_samples_contacts_table.write_html(f'{args.output_dir}/plots/samples_contacts_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        contacts_message = '''Table of mutated residues with other residue interactions, and their interaction type e.g. H-bond, salt-bridge, "contact".
                            A full table of interacting residue mutations in each sample can be found <a href="plots/samples_contacts_table.html">here</a> and in <code>spear_annotation_summary.tsv</code>'''

    else:
        contacts_table_plt = "No mutations in contacting residues present in any samples."
        contacts_message = ""

    #MAKING A SCORES TABLE COLOURED WHERE SAMPLE SCORE SUM EXCEEDS BASELINE
    scores_cols = baseline_scores.columns.tolist()
    non_displayed_scores = ['total_variants', 'total_residue_variants','consequence_type_variants', 'region_residues', 'domain_residues', 'feature_residues', 'ACE2_contact_counts', 'ACE2_contact_score', 'trimer_contact_counts','trimer_contact_score', 'barnes_class_variants','bloom_ACE2_wuhan_max', "bloom_ACE2_wuhan_min", 'bloom_ACE2_BA1_min', 'bloom_ACE2_BA1_max', 'bloom_ACE2_BA1_min', 'bloom_ACE2_BA2_max', 'bloom_ACE2_BA2_min', 'VDS_max', 'VDS_min', 'serum_escape_max', 'serum_escape_min', 'cm_mAb_escape_all_classes_max','cm_mAb_escape_all_classes_min','mAb_escape_all_classes_max', 'mAb_escape_all_classes_min', 'mAb_escape_class_1_max', 'mAb_escape_class_1_min', 'mAb_escape_class_2_max', 'mAb_escape_class_2_min', 'mAb_escape_class_3_max', 'mAb_escape_class_3_min', 'mAb_escape_class_4_max', 'mAb_escape_class_4_min', 'BEC_RES_max', 'BEC_RES_min', 'BEC_RES_sum']
    displayed_scores_cols = [score for score in scores_cols if score not in non_displayed_scores]
    sample_scores = scores_summary[displayed_scores_cols]
    sample_scores = sample_scores.replace("", np.nan).dropna(axis=1, how = "all") #remove empty cols from table to be displayed (do this later for the graph table to allow subtraction of baseline array)
    displayed_scores_cols = [score for score in displayed_scores_cols if score in sample_scores.columns.tolist()]
    actual_scores_cols = [score for score in displayed_scores_cols if score != "sample_id"]
    if sample_scores[[score for score in actual_scores_cols if score not in ("displayed_dropout","lineage")]].isna().all().all() == False:
        if "cm_mAb_escape_all_classes_sum" in sample_scores.columns:
            sort_col = "cm_mAb_escape_all_classes_sum"
        else:
            sort_col = "sample_id"
        sample_scores = sample_scores.replace("", np.nan).sort_values(by = sort_col, ascending = False)
        sample_scores = pd.concat([baseline_scores.loc[baseline_scores["sample_id"] == args.baseline, displayed_scores_cols], sample_scores])
        sample_scores = sample_scores.reset_index(drop = True)
        #if baseline is also in sample set this can cause problems for the subtraction of a single numpy array from dataframe. Baseline will always be first index value so can fallback on this if multiple matches
        if len(sample_scores.loc[sample_scores["sample_id"] == args.baseline, "sample_id"]) > 1:
            baseline_relative_sample_scores = sample_scores[actual_scores_cols] - sample_scores.loc[0, actual_scores_cols].fillna(0).values.squeeze()
        else:    
            baseline_relative_sample_scores = sample_scores[actual_scores_cols] - sample_scores.loc[sample_scores["sample_id"] == args.baseline, actual_scores_cols].fillna(0).values.squeeze()
        
        baseline_relative_sample_scores = baseline_relative_sample_scores.round(4)
        
        sample_scores = sample_scores[sample_scores.index.isin(baseline_relative_sample_scores.index)].reset_index(drop = True)
        baseline_relative_sample_truths = np.greater(baseline_relative_sample_scores.fillna(-1).to_numpy(), 0) #all others should be coloured if greater than baseline
        baseline_relative_sample_colours = np.where(baseline_relative_sample_truths == True, "rgb(250,180,174)", "rgb(179,205,227)")
        baseline_relative_sample_colours[0,:] = np.array(['lavender'] * np.shape(baseline_relative_sample_colours)[1])

        sample_id_colours = np.array([['lavender'] * np.shape(baseline_relative_sample_colours)[0]]).T
        baseline_relative_sample_colours = np.append(sample_id_colours, baseline_relative_sample_colours, axis = 1)
        baseline_relative_sample_colours_df = pd.DataFrame(baseline_relative_sample_colours, columns = displayed_scores_cols)
        
        baseline_relative_sample_colours_df.set_index(sample_scores["sample_id"], drop = False, inplace = True)
        sample_scores[actual_scores_cols] = sample_scores[actual_scores_cols].round(2).fillna("").astype("str")
        
        labels = {
            "sample_id" : "Sample ID", 
            "lineage" : "Pango-Lineage",
            "VDS_mean" : "Vibrational Difference Score",
            "bloom_ACE2_wuhan_mean" : "Bloom ACE2 (Wuhan)",
            "bloom_ACE2_BA1_mean" : "Bloom ACE2 (BA.1)",
            "bloom_ACE2_BA2_mean" : "Bloom ACE2 (BA.2)",
            "serum_escape_sum" : "Serum Escape",
            "mAb_escape_all_classes_sum" : "mAb Escape",
            "cm_mAb_escape_all_classes_sum" : "Class Masked mAb Escape",
            "mAb_escape_class_1_sum" : "mAb Escape Class 1",
            "mAb_escape_class_2_sum": "mAb Escape Class 2",
            "mAb_escape_class_3_sum": "mAb Escape Class 3",
            "mAb_escape_class_4_sum": "mAb Escape Class 4",
            "BEC_EF_sample" : "BEC Escape Factor",
            "displayed_dropout" : "Quality Warnings"}

        #now add in the S gene dropout and Global N percentage columns 
        #will have to add a colour column to the np color array too baseline_relative_sample_colours_df, add a column to displayed_scores_cols too
        if os.path.basename(args.n_perc) == "spear_score_summary.tsv":
            sample_scores["displayed_dropout"] = ""
        else:
            n_info = pd.read_csv(args.n_perc)
            n_info.loc[n_info["s_n_contig"] >= s_contig, "s_n_contig_display"] = "!"
            n_info.loc[n_info["rbd_n"] >= rbd_n, "rbd_n_display"] = "^"
            n_info.loc[n_info["s_n"] >= s_n, "s_n_display"] = "#"
            n_info.loc[n_info["global_n"] >= global_n, "global_n_display"] = "*"
            n_info[["s_n_contig_display", "s_n_display", "global_n_display", "rbd_n_display"]] = n_info[["s_n_contig_display", "s_n_display", "global_n_display", "rbd_n_display"]].fillna("")
            n_info["displayed_dropout"] = n_info["s_n_contig_display"] + n_info["rbd_n_display"] + n_info["s_n_display"] + n_info["global_n_display"]
            sample_scores_baseline = sample_scores.iloc[0].copy()
            sample_scores_baseline["displayed_dropout"] = ""
            sample_scores_samples = sample_scores.iloc[1:].copy()
            sample_scores_samples = pd.merge(left = sample_scores_samples, right = n_info[["sample_id", "displayed_dropout"]], on = "sample_id", how = "left").fillna("")
            sample_scores = pd.concat([sample_scores_baseline.to_frame().T, sample_scores_samples])
        sample_scores.set_index("sample_id", drop = False, inplace = True)
        sample_scores.index.name = "index"

        n_info_colours = np.where(sample_scores["displayed_dropout"] != "", "rgb(250,180,174)", "rgb(179,205,227)")
        n_info_colours[0] = "lavender"
        lineages_sub = lineages[["taxon", "lineage"]].rename(columns = {"taxon" : "sample_id"})
        sample_scores = pd.merge(sample_scores, lineages_sub, on = "sample_id" , how = "left")
        sample_scores.index = sample_scores["sample_id"]
        sample_scores = sample_scores.rename_axis("index")

        displayed_scores_cols.append("displayed_dropout")
        baseline_relative_sample_colours_df["displayed_dropout"] = n_info_colours
        displayed_scores_cols.insert(1,"lineage")
        baseline_relative_sample_colours_df["lineage"] = "lavender"
        baseline_relative_sample_colours = np.concatenate([baseline_relative_sample_colours, np.reshape(n_info_colours, (-1,1))], axis = 1)
        scores_table = go.Figure(data=[go.Table(
            header=dict(values= [labels[col] for col in displayed_scores_cols],
                        fill_color='paleturquoise',
                        align='center'),
            cells=dict(values= [sample_scores[x] for x in displayed_scores_cols],
                    fill_color=[baseline_relative_sample_colours_df[x] for x in displayed_scores_cols],
                    align='center'))
                    ])  
        buttons = []
        for score in displayed_scores_cols:
            if score == "sample_id" or score == "lineage":
                asc = True
            else:
                asc = False
            baseline_scores = sample_scores.iloc[[0]]
            baseline_colours = baseline_relative_sample_colours_df.iloc[[0]]
            samples_scores = sample_scores.iloc[1:, :]
            samples_colours = baseline_relative_sample_colours_df.iloc[1:, :]
            samples_scores = samples_scores.replace("", np.nan).sort_values(by = score, ascending = asc ).replace(np.nan, "")
            samples_colours = samples_colours.reindex(samples_scores.index)
            
            sorted_scores = pd.concat([baseline_scores, samples_scores])
            sorted_colours = pd.concat([baseline_colours, samples_colours])

            buttons.append(dict(
                    label = labels[score],
                    method = 'restyle',
                    args = [
                        {"cells": {
                            "values": [sorted_scores[x] for x in displayed_scores_cols], 
                            "fill": dict(color = [sorted_colours[x] for x in displayed_scores_cols])}}]))

        scores_table.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = displayed_scores_cols.index(sort_col),
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

        scores_table.update_layout({"paper_bgcolor":'rgba(0,0,0,0)', "margin" : dict(r=5, l=5, t=5, b=5), "autosize" : True})
        scores_table_plt = offline.plot(scores_table, output_type='div', include_plotlyjs = False , config = {'displaylogo': False})
        scores_table.write_html(f'{args.output_dir}/plots/scores_table.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        score_table_message = '''For a full screen view of this table see <a href="plots/scores_table.html">here</a>. Source data used to produce this table can be found in the file <code>spear_score_summary.tsv</code>'''

        table = Table(show_header=True, header_style="bold magenta", title = "Per Sample Scores Summary", caption = "Quality warnings: ! - Spike N contig (default 150nt)  ;  ^ - Spike RBD N content (default 12nt)  ;  * - Global N percentage (default > half N percentage cutoff) ;  # - Spike N percentage (default > 5%)", caption_justify = "center")
        pango_lineage = sample_scores.pop("lineage")
        sample_scores.insert(1, "lineage", pango_lineage)
        pango_lineage_colours = baseline_relative_sample_colours_df.pop("lineage")
        baseline_relative_sample_colours_df.insert(1, "lineage", pango_lineage_colours)
        for column in sample_scores.columns:
            table.add_column(labels[column])
        cli_baseline_relative_sample_truths = np.where(np.isin(baseline_relative_sample_colours_df,["rgb(179,205,227)", "lavender"]), False, True)
        for x, y in zip(sample_scores.values, cli_baseline_relative_sample_truths):
            row_value = []
            for (a, b) in zip(x,y):
                colour = "#FF0000" if b else "#2CBDC9"
                value = a if a == a else ""
                bold = "bold" if b else "default"
                row_value.append((f'[{bold} {colour}]{value}'))
            table.add_row(*row_value)
    else:
        table = "No variants were detected in scoring regions of Spike"
        scores_table_plt = "No variants were detected in scoring regions of Spike"
        score_table_message = ""

    #MAKING THE INTERACTIVE PLOTS: 
    scores_cols = ["bloom_ACE2_wuhan", "bloom_ACE2_BA1", "bloom_ACE2_BA2", "VDS","serum_escape", "mAb_escape_all_classes", "cm_mAb_escape_all_classes", "mAb_escape_class_1", "mAb_escape_class_2", "mAb_escape_class_3", "mAb_escape_class_4", "BEC_RES"]
    scores_z_max = {"bloom_ACE2_wuhan" : 4.84, "bloom_ACE2_BA1" : 4.84, "bloom_ACE2_BA2" : 4.84, "VDS": 0.712636025 , "serum_escape" : 1 , "mAb_escape_all_classes" : 1, "cm_mAb_escape_all_classes" : 1, "mAb_escape_class_1" : 1, "mAb_escape_class_2" : 1, "mAb_escape_class_3" : 1, "mAb_escape_class_4" : 1, "BEC_RES" : 1}
    scores_z_min = {"bloom_ACE2_wuhan" : -4.84, "bloom_ACE2_BA1" : -4.84, "bloom_ACE2_BA2" : -4.84, "VDS" : -0.712636025 ,"serum_escape" : 0 , "mAb_escape_all_classes" : 0, "cm_mAb_escape_all_classes" : 0, "mAb_escape_class_1" : 0, "mAb_escape_class_2" : 0, "mAb_escape_class_3" : 0, "mAb_escape_class_4" : 0, "BEC_RES" : 0}
    scores_z_mid = {"bloom_ACE2_wuhan" : 0, "bloom_ACE2_BA1" : 0, "bloom_ACE2_BA2" : 0, "VDS" : 0, "serum_escape" : 0.5, "mAb_escape_all_classes" : 0.5, "cm_mAb_escape_all_classes" : 0.5, "mAb_escape_class_1" : 0.5, "mAb_escape_class_2" : 0.5, "mAb_escape_class_3" : 0.5, "mAb_escape_class_4" : 0.5, "BEC_RES" : 0.5}
    scores_title = {"bloom_ACE2_wuhan" : "Bloom ACE2 (Wuhan)", "bloom_ACE2_BA1" : "Bloom ACE2 (BA1)", "bloom_ACE2_BA2" : "Bloom ACE2 (BA2)", "VDS" : "Vibrational Difference Score","serum_escape" : "Serum Escape", "mAb_escape_all_classes" : "mAb Escape", "cm_mAb_escape_all_classes" : "Class Masked mAb Escape", "mAb_escape_class_1" : "mAb Escape Class 1", "mAb_escape_class_2": "mAb Escape Class 2", "mAb_escape_class_3": "mAb Escape Class 3", "mAb_escape_class_4": "mAb Escape Class 4", "BEC_RES" : "BEC Residue Escape Score "}
    scores_color_scales = {"bloom_ACE2_wuhan" : "plasma", "bloom_ACE2_BA1" : "plasma", "bloom_ACE2_BA2" : "plasma", "VDS" : "rdbu","serum_escape" : "hot_r", "mAb_escape_all_classes" : "hot_r", "cm_mAb_escape_all_classes" : "hot_r", "mAb_escape_class_1" : "hot_r", "mAb_escape_class_2": "hot_r", "mAb_escape_class_3": "hot_r", "mAb_escape_class_4": "hot_r", "BEC_RES" : "purd_r"}
    
    respos_df = pd.read_csv(f'{args.data_dir}/product_mapping.csv')
    orf_boxes = []
    orf_labels = []
    for product in products:
        orf_boxes.append(add_orf_shape(respos_df.loc[respos_df["product"] == product, "position"].min(), respos_df.loc[respos_df["product"] == product, "position"].max(), product, orf1ab, structural, accessory,product_colours))
        orf_labels.append(add_orf_label(respos_df.loc[respos_df["product"] == product, "position"].min(), respos_df.loc[respos_df["product"] == product, "position"].max(), product, orf1ab, structural, accessory))

    count = 0
    if args.product_plots == True:
        table_html_string = []
        table_header = '''<table class="table table-responsive text-center table-hover"><tr><td>Sample ID</td><td>Product Plot</td><td>AA Mutations</td></tr>''' #style="width: auto;"
        table_html_string.append(table_header)
        for sample in track(annotation_summary["sample_id"].unique().tolist(), description="Saving per sample mutation plots..."):
            sample_annotation_summary = annotation_summary.loc[annotation_summary["sample_id"] == sample].copy()
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
            for score in scores_cols:
                if score == "cm_mAb_escape_all_classes":
                    visible = True
                else:
                    visible = False
                sample_orf_plot.add_trace(go.Scatter(
                    {
                        "x" : merged_sample_occurences["position"], 
                        "y" : merged_sample_occurences["axis_pos"], 
                        "text" : merged_sample_occurences["residues"],
                        "mode" : "markers+text", 
                        'visible' : visible,
                        "textposition" :  merged_sample_occurences["textpos"],
                        "hovertemplate" : merged_sample_occurences["residues"] + "<br>" + merged_sample_occurences[score].fillna("").astype("str") + "</br><extra></extra>",
                        "marker" : {
                            "size" : 12,
                            "color" : merged_sample_occurences[score],
                            "colorbar" : dict(title = ""),
                            "colorscale" : scores_color_scales[score],
                            "cmid": scores_z_mid[score],
                            "cmax": scores_z_max[score],
                            "cmin": scores_z_min[score],
                            "line" : {"width" : 2 , "color": 'DarkSlateGrey'}
                        }
                    }))

            buttons = []
            trace_list = [True] * len(scores_cols)
            count = 0
            for score in scores_cols:
                trace_list = [False] * len(scores_cols)
                trace_list[count] = True
                buttons.append(dict(
                        label = scores_title[score],
                        method = 'update',
                        args = [
                            {'visible': trace_list}]))
                count += 1
            sample_orf_plot.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',xaxis=dict(showgrid=False, showline = False, zeroline = False),yaxis=dict(showgrid=False, showline = False, zeroline = False)) #title_text='WHAT SHOULD THIS TITLE BE ? ', title_x=0.5, 
            sample_orf_plot.update_layout(
                annotations = orf_labels,
                shapes = combo, 
                title= {"text" : f'{sample} Product Plot', "x" : 0.5},
                xaxis = {"range" : [0,10000], "fixedrange" : True, "showticklabels" : False, "showgrid" : False},
                yaxis = {"range" : [-14,18], "fixedrange" : True, "showticklabels" : False, "showgrid" : False})
            
            # Add dropdown
            sample_orf_plot.update_layout(
                updatemenus=[
                    dict(
                        buttons=list([
                            dict(
                                args=["mode", "markers+text"],
                                label="Marker + mutation",
                                method="restyle"),
                            dict(
                                args=["mode", "markers"],
                                label="Marker only",
                                method="restyle"
                            ),
                            dict(
                                args=["mode", "text"],
                                label="Mutation only",
                                method="restyle"
                            ),
                        ]),
                        active = 0,
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.16,
                        xanchor="left",
                        y=1.075,
                        yanchor="top"
                    ),
                    dict(
                        buttons=list([
                            dict(
                                args=[{"shapes" : combo}],
                                label="Show lines",
                                method="relayout"),
                            dict(
                                args=[{"shapes" : orf_boxes}],
                                label="Hide lines",
                                method="relayout"),
                        ]),
                        active = 0,
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.01,
                        xanchor="left",
                        y=1.075,
                        yanchor="top"
                    ),
                    dict(
                        buttons= buttons,
                        active = scores_cols.index("cm_mAb_escape_all_classes"),
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.16,
                        xanchor="left",
                        y=1.01,
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
                                    {"xaxis" : {"range" :[7000,8400], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[3.3,16], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="Spike Protein",
                                method="relayout"),
                            dict(
                                args = [
                                    {"xaxis" : {"range" :[8658,9719], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[3.3,16], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="Structural (non-Spike)",
                                method="relayout"),
                            dict(
                                args = [
                                    {"xaxis" : {"range" :[8383,9757], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[-14,-2], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="Accessory proteins",
                                method="relayout"),
                            dict(
                                args = [
                                    {"xaxis" : {"range" :[0,4393], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[-12,12], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}},
                                    {"annotations" : {"font" : {"size" : 50}}}
                                ],
                                label="pp1a",
                                method="relayout"),
                            dict(
                                args = [
                                    {"xaxis" : {"range" :[4394,7109], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}, "yaxis" : {"range" :[-12,12], "fixedrange" : True, "showticklabels" : False, "showgrid" : False}}
                                ],
                                label="pp1ab only",
                                method="relayout"),
                        ]),
                        active = 0,
                        bgcolor = "white",
                        direction="down",
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.01,
                        xanchor="left",
                        y=1.01,
                        yanchor="top"
                    ),
                ]
            )

            sample_orf_plot.write_html(f'{args.output_dir}/plots/product_plots/{sample}_product_plot.html', include_plotlyjs = f'../plotly/plotly-2.8.3.min.js')
            sample_total_variants = merged_sample_occurences["residues"].count()
            sample_table_row = f'<tr><td>{sample}</td><td><a href="plots/product_plots/{sample}_product_plot.html">{sample} Product Plot</a></td><td>{sample_total_variants}</td>'
            table_html_string.append(sample_table_row)
        table_html_string.append("</table>")
        table_html_string = "".join(table_html_string)
        table_card = '''
                <div class="row justify-content-center mt-2">
                    <div class="col-auto">
                        <div class="card text-center">
                            <div class="card-header" id="headingOne">
                                <h5 class="mb-2">
                                    <button class="btn btn-link" data-toggle="collapse" data-target="#collapseOne" aria-expanded="false" aria-controls="collapseOne">
                                    ORF Product Plots
                                    </button>
                                </h5>
                            </div>
                            <div id="collapseOne" class="collapse hide" aria-labelledby="headingOne">
                                <div class="card-body">
                                    <div class="table-responsive">
                                    ''' + table_html_string + '''
                                    </div>
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
    anno_merge = pd.merge(anno_merge, lineages[["taxon", "lineage"]], left_on = "sample_id", right_on = "taxon", how = "left")
    anno_merge.set_index("residues")
    if anno_merge.loc[anno_merge["lineage"].isna() == False, "lineage"].isin(["", np.nan]).all():
        anno_merge["text_var"] = anno_merge["sample_id"] + ": " + anno_merge["residues"]
        anno_merge["lineage_sample"] = anno_merge["sample_id"]
    else:
        anno_merge["text_var"] = anno_merge["sample_id"] + "(" + anno_merge["lineage"] + "): " + anno_merge["residues"]
        anno_merge["lineage_sample"] = anno_merge["lineage"] + anno_merge["sample_id"]

    lineage_categories = anno_merge.loc[anno_merge["lineage_sample"].isna() == False, "lineage_sample"].to_list()
    lineage_categories = sorted(lineage_categories)
    displayed_scores = []
    if anno_merge[scores_cols].isna().all().all() == False:
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
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 14) & (anno_merge["respos"] <= 913),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : heatmap_text,
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "rdbu",
                            "name" : "VDS",
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                        }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 14) & (anno_merge["residue-position"] <= 913),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : heatmap_all_text,
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "rdbu",
                            "name" : "VDS",
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                        }))
                elif score == "bloom_ACE2_wuhan":
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma", 
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False, 
                        }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma",
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                        }))
                elif score == "bloom_ACE2_BA1":
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma", 
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False, 
                        }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma",
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                        }))
                elif score == "bloom_ACE2_BA2":
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma", 
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False, 
                        }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            "colorscale" : "plasma",
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                        }))

                elif score == "cm_mAb_escape_all_classes":                
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "hot_r",
                            "name" : "cm_mAb_escape_all_classes", 
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                            }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "hot_r",
                            "name" : "cm_mAb_escape_all_classes",
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                            }))
                elif score == "BEC_RES":                
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "name" : "BEC_RES",
                            'colorscale' : "purd_r",
                            }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "name" : "BEC_RES",
                            'colorscale' : "purd_r",
                            }))
                else:
                    heatmap.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"respos"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["respos"] >= 331) & (anno_merge["respos"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "hot_r",
                            "name" : score, 
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                            }))
                    heatmap_all.add_trace(go.Heatmap(
                        {
                            'z': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531), score].values.tolist(),
                            'x': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"lineage_sample"].astype("category").values.tolist(),
                            'y': anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"residue-position"].astype("Int64").astype("str").values.tolist(),
                            'text' : anno_merge.loc[(anno_merge["product"] == "surface glycoprotein") & (anno_merge["residue-position"] >= 331) & (anno_merge["residue-position"] <= 531),"text_var"].values.tolist(),
                            'texttemplate' : "%{text}",
                            'hovertemplate' : '%{text} <br>Score: %{z}<extra></extra>',
                            'visible' : False,
                            "colorscale" : "hot_r", 
                            "zmin" : scores_z_min[score],
                            "zmax" : scores_z_max[score],
                            "zmid" : scores_z_mid[score]
                            }))

        if "cm_mAb_escape_all_classes" not in displayed_scores:
            active = 0 #for buttons set to first button
            heatmap.update_traces(visible = True,
                        selector=dict(name=displayed_scores[0])) #set to first non "sample_id" score col i.e. 1  
            heatmap_all.update_traces(visible = True,
                        selector=dict(name=displayed_scores[0])) #set to first non "sample_id" score col i.e. 1  
        else:
            active = displayed_scores.index("cm_mAb_escape_all_classes")
            heatmap.update_traces(visible = True,
                        selector=dict(name="cm_mAb_escape_all_classes"))
            heatmap_all.update_traces(visible = True,
                        selector=dict(name="cm_mAb_escape_all_classes"))
        layout = {"title" : dict(text = "Class Masked mAb Escape", x = 0.5), 
            "xaxis" : {"title": "Sample" , "showticklabels" : False, "showgrid" : False},
            "yaxis" : {"title": "Residue Position" ,"tickformat": '.0f', "showgrid" : False}}
        heatmap.update_layout(layout)
        heatmap_all.update_layout(layout)
        heatmap.update_xaxes(type="category")
        heatmap.update_xaxes(categoryorder='array', categoryarray= lineage_categories)
        #heatmap_all.update_yaxes(categoryorder="category ascending")
        
        trace_list = [True] * len(displayed_scores)
        count = 0
        buttons = []

        for score in displayed_scores:
            trace_list = [False] * len(displayed_scores)
            trace_list[count] = True
            buttons.append(dict(
                            label = scores_title[score],
                            method = 'update',
                            args = [
                                {'visible': trace_list},
                                {'title': dict(text = scores_title[score], x = 0.5)},
                                {'yaxes_type' : 'category'}]))
            count += 1
        heatmap.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor = "rgb(246,246,246)")  
        heatmap_all.update_layout(paper_bgcolor = 'rgba(0,0,0,0)', plot_bgcolor = "rgb(246,246,246)")
        heatmap.update_layout(
            updatemenus=[
                dict(
                    buttons=buttons,
                    active = active,
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
                    active = active,
                    bgcolor = "white",
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x = 0.01,
                    y = 1.3,
                    xanchor="left",
                    yanchor="top")
                    ])

        heatmap_html = offline.plot(heatmap,output_type='div', include_plotlyjs = False, config = {'displaylogo': False})

        heatmap.write_html(f'{args.output_dir}/plots/mutated_residues_heatmap.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        heatmap_all.write_html(f'{args.output_dir}/plots/all_residues_heatmap.html', include_plotlyjs=f'plotly/plotly-2.8.3.min.js')
        heatmap_message = f'For a full screen view of the current plot see <a href="plots/mutated_residues_heatmap.html">here</a>, and for a fullscreen heatmap across all residues (easier comparison between reports), see <a href="plots/all_residues_heatmap.html">here</a>.'
    else:
        heatmap_html = "<p>No variants were detected in scoring regions of Spike</p>"
        heatmap_message  = ""
    
    #################### HTML FORMATTING #####################
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
            <style>body{ margin:0 100; background:white; }</style>
        </head>
        <body>
            <nav class="navbar navbar-light bg-light navbar-expand-md">
                <div class="navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2">
                    <ul class="navbar-nav mr-auto">
                        <li class="nav-item active">
                            <a class="nav-link" href="#"><img src="images/SPEAR_smallest.png" class="img-fluid" alt="Responsive image"></a>
                        </li>
                    </ul>
                </div>
                <div class="mx-auto order-0">
                    <a class="navbar-brand text-center" href="#"><span class="mb-0 h1">Summary Report</span></a>
                    <button class="navbar-toggler" type="button" data-toggle="collapse" data-target=".dual-collapse2">
                        <span class="navbar-toggler-icon"></span>
                    </button>
                </div>
                <div class="navbar-collapse collapse w-100 order-3 dual-collapse2">
                    <ul class="navbar-nav ml-auto">
                    </ul>
                </div>
            </nav>
            <div class="container-fluid">
                <div class = "row mt-2">
                    <div class = "col-12">
                        <div class="card">
                            <div class="card-body">
                            <p class="p1 text-justify">From a total of ''' + str(spear_qc_info["input_samples"].values[0]) + ''' input samples, ''' + str(qc_samples) + ''' passed QC and were annotated by SPEAR. A total of ''' + str(total_genomic_variants) + ''' nucleotide variants were identified across all samples, which resulted in ''' + str(total_residue_variants) + ''' AA missense, deletion or insertions which are listed and evaluated below. </p>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="row mt-2">
                    <div class="col-12">
                        <div class="card">
                            <div class="card-header">Amino Acid Variants</div>
                            <div class="card-body">
                                <div>
                                ''' + residues_table_plt + '''
                                </div>
                            </div>
                            <div class="card-footer">Amino acid missense, deletion, or insertion mutations across all samples. Sorted by count of mutation, use dropdown to sort by genomic position. Nucleotide level variants may be repeated where there are multiple AA changes associated with a single genomic event. A table summarising nucleotide level variants can be found <a href="plots/nt_variants_table.html">here</a>.</div>
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
                                Summary heatmap of all residue positions for which there is at least one mutation across all samples. 
                                Individual sample IDs and mutations are plotted on datapoint, for large samples sets these become visible at higher zoom levels. Residues are labelled in the following format: sample_id(pango_lineage):mutation.   
                                Hover text will show z: selected score for sample-residue and sample ID and mutation.  
                                For a description of these scores see <a href="https://github.com/m-crown/SPEAR#scores">Table 3</a> in the SPEAR README.  
                                ''' + heatmap_message + '''
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
                            <div class = "card-footer">
                            Summarised scores per sample (sum across sample), cells with values higher than selected baseline are highlighted. 
                            Selected baseline is always shown in the top row, table is sorted by the cm mAb escape all classes sum column by default, 
                            the drop down can be used to sort on other scores. For a description of these scores see <a href="https://github.com/m-crown/SPEAR/blob/main/docs/Table4.md#spear-score-summary">Table 4</a> in the SPEAR README. ''' + score_table_message + '''<br>
                            Quality warnings: ! - Spike N contig (default 150nt)  ;  ^ - Spike RBD N content (default 12nt)  ;  * - Global N percentage (default > half N percentage cutoff) ;  # - Spike N percentage (default > 5%)</br>
                            </div>
                        </div>
                    </div> 
                </div>
                <div class = "row mt-2">
                    <div class="col-6">
                        <div class="card">
                            <div class="card-header">
                                Protein Domain Mutations
                            </div>
                            <div class="card-body">
                                ''' + domain_table_plt + '''
                            </div>
                            <div class = "card-footer">
                            ''' + domain_message + '''
                            </div>
                        </div>
                    </div>
                    <div class="col-6"> 
                        <div class="card">
                            <div class="card-header">
                                Protein Feature Mutations
                            </div>
                            <div class="card-body">
                                ''' + feature_table_plt + '''
                            </div>
                            <div class = "card-footer">
                            ''' + feature_message + '''
                            </div>
                        </div> 
                    </div>    
                </div>
                <div class = "row mt-2">
                    <div class="col-6">
                        <div class="card">
                            <div class = "card-header">
                                Mutated Interacting Residues
                            </div>
                            <div class = "card-body">
                                ''' + contacts_table_plt + '''
                            </div>
                            <div class = "card-footer">
                            ''' + contacts_message + '''
                            </div>
                        </div>
                    </div>
                    <div class="col-6">
                        <div class="card">
                            <div class = "card-header">
                                Mutated Furin Site Residues
                            </div>
                            <div class = "card-body">
                                ''' + furin_table_plt + '''
                            </div>
                            <div class = "card-footer">
                            ''' + furin_message + '''
                            </div>
                        </div>
                    </div>
                </div>
                <div class = "row mt-2">
                    <div class = "col-6">
                        <div class="card">
                            <div class = "card-header">
                                Mutated NTD Residues
                            </div>
                            <div class = "card-body">
                                ''' + ntd_table_plt + '''
                            </div>
                            <div class = "card-footer">
                            ''' + ntd_message + '''
                            </div>
                        </div>
                    </div>
                </div>
                ''' + table_card + '''    

            <footer class="page-footer font-small teal pt-4 border mt-2 mb-2">
                <div class="container-fluid">
                    <div class="row">
                        <div class="col-4 text-left">
                            <p> Generated On: ''' + report_date + '''</p>
                        </div>
                        <div class="col-4 text-center">
                            <p> ''' + spear_version + ''' </p>
                        </div>
                        <div class="col-4 text-right">
                            <p></p>
                        </div>
                    </div>                   
                    <div class="row">
                        <div class="col-md-12 mt-md-0 mt-3 text-center">
                            <p><a href="https://github.com/m-crown/SPEAR">SPEAR</a> was developed by Matt Crown and Matt Bashton at <a href="https://www.northumbria.ac.uk/">Northumbria University UK</a> from work funded by <a href="https://www.cogconsortium.uk/">COG-UK</a>. Please post bugs, issue and feature requests on <a href="https://github.com/m-crown/SPEAR/issues">GitHub</a>. </p>
                        </div>
                    </div>
                    <div class="row">
                        <div class="col">
                            <h6>''' + spear_params +'''</h6>
                        </div>
                    </div>
                    <div class="row">
                        <div class="col">
                            <h6>''' + pangolin_params +'''</h6>
                        </div>
                    </div>
                    <div class="row">
                        <div class="col">
                            <h6>Pangolin versions: ''' + pangolin_versions +'''</h6>
                        </div>
                    </div>
                </div>
            </footer>
        </div>
        <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
        <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
        </body>
    </html>'''

    # <div class="col-6">
    #     <div class="card">
    #         <div class="card-header">Nucleotide Variants</div>
    #         <div class="card-body">
    #             <div>
    #             ''' + variants_table_plt + '''
    #             </div>
    #         </div>
    #         <div class="card-footer">Genomic variants detected across all samples. Sorted by variant count, use dropdown to sort by genomic position.</div>
    #     </div>
    # </div>

    copyfile(f'{args.images_dir}/SPEAR_smallest.png', f'{args.output_dir}/images/SPEAR_smallest.png')
    f = open(f'{args.output_dir}/report.html','w')
    f.write(html_string)
    f.close()
    console.print(table)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--product_plots', default=False, action='store_true',
        help = "Output per sample ORF plots.")
    parser.add_argument('--n_perc', metavar='n_perc.csv', type=str,
        help='Filename for %N and S gene dropout')
    parser.add_argument('score_summary', metavar='spear_score_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF')
    parser.add_argument('annotation_summary', metavar='spear_annotation_summary.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF')
    parser.add_argument('baseline_scores', metavar='baseline.tsv', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF')
    parser.add_argument('images_dir', metavar='spear_images/', type=str,
        help='Directory for spear images to be copied from') 
    parser.add_argument('scripts_dir', metavar='$CONDA_PREFIX/bin', type=str,
        help='Directory for spear scripts')
    parser.add_argument('data_dir', metavar='$CONDA_PREFIX/data', type=str,
        help='Directory for spear data')        
    parser.add_argument('output_dir', metavar='report/', type=str,
        help='Filename for SnpEff summarised and residue annotated VCF')
    parser.add_argument('baseline', metavar='Omicron', type=str,
        help='lineage for baseline')
    parser.add_argument('pangolin_report', metavar='lineage_report.csv', type=str,
        help='pangolin lineage report for organising samples')
    parser.add_argument('pangolin_command', metavar='pangolin_command.txt', type=str,
        help='pangolin command which was run')
    parser.add_argument('spear_params', metavar='spear_params.csv', type=str,
        help='spear params')
    parser.add_argument('spear_qc_info', metavar='spear_qc_info.tsv', type=str,
        help='spear qc info')
    args = parser.parse_args()

    summary_report(args = args)

if __name__ == "__main__":
    main()
