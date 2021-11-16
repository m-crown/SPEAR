#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
    File name: go_propa.py
    Author: Matthew Crown
    Date created: 22/10/2021
    Python Version: 3.8.3

    Tool to convert the absolute gene product coordinate system of gofasta sam 
    variants to a mature peptide relative position for specified polyproteins 
    e.g. for Sars-CoV-2 ORF1ab and ORF1ab polyproteins.
'''
# ==============================================================================
from Bio import SeqIO
import csv
import copy
import argparse
import json

def variants_parser(variants_list):
    '''
    Takes the ouput CSV from gofasta sam variants and stores samples in a nested
    dictionary containing lists of mutations for each gene.
    '''
    samples = {}
    with open(variants_list, 'r') as file:
        csvreader = csv.reader(file)
        header = next(csvreader)
        for sample in csvreader:
            if sample[1] != '':
                samples[sample[0]] = {}
                #split out pipe delimited variants
                for var in sample[1].split("|"):
                    #split variant call into gene and position 
                    gene, pos = var.split(":") 
                    #check if the gene already exists in nested dictionary keys
                    if gene in samples[sample[0]].keys(): 
                        #deduplicate the output from gofasta sam variants - particularly needed for orf1ab variants which are in orf1a and orf1ab
                        if pos not in samples[sample[0]][gene].keys(): 
                            #append position to exiting list for this gene
                            samples[sample[0]][gene][pos] = {} 
                    else:
                        #if gene doesnt already exist in nested dictionary create the key and set values to be a dict including the position
                        samples[sample[0]][gene] = {pos: {}} 
    return samples

def get_variant_info(samples, coords):
    '''
    Function takes a list of variants with polyprotein positions and converts 
    them to matured peptide positions; also returns name of the matured peptide 
    the variant is in
    '''
    for sample in samples.values():        
        for gene_name, variants_list in sample.items():
            if gene_name != "synSNP":
                for variant_name, variant_info in variants_list.items():
                    variant_info["product"] = coords[gene_name]["product"]
                    if "multi-product" in coords[gene_name]:
                        placeholder_names = []
                        placeholder_pos = []
                        for k,v in coords[gene_name]["multi-product"].items():
                            if v["sub-position"][0] <= int(variant_name[1:-1]) -1 <= v["sub-position"][1]:
                                #using a placeholder for the edge case situation where multiproducts overlap and variant is in the overlap. First collect all of the matching names/pos into placeholder
                                placeholder_names.append(k) 
                                placeholder_pos.append(f'{variant_name[0]}{int(variant_name[1:-1]) - v["sub-position"][0]}{variant_name[-1]}')
                        #If there are more than one item in the placeholder names list join together names and positions with /, else set name/pos to the only matching gene/product.
                        if len(placeholder_names) > 1:
                            name = "/".join(placeholder_names)
                            if len(set(placeholder_pos)) == 1:
                                pos = placeholder_pos[0]
                        else:
                            name = placeholder_names[0]
                            pos = placeholder_pos[0]
                        #set the variant info multi-product name and subposition.
                        variant_info["multi-product"] = name
                        variant_info["sub-position"] = pos
    return samples

def flatten_variants(samples):
    flattened_vars = []
    flattened_vars.append(["sample","gene", "position", "mat-position", "product"])
    for sample , genes in samples.items():
        for gene, desc in genes.items():
            for var, description in desc.items():
                if "multi-product" in description:
                    variant_description = [sample, gene, var, description["sub-position"], description["multi-product"]]
                    flattened_vars.append(variant_description)
                elif gene == "synSNP" or gene == "del" or gene == "ins": #the three possible outputs from gofasta which are genomic coords rather than peptide
                    variant_description = [sample,gene,var,"",""]
                    flattened_vars.append(variant_description)
                else:
                    variant_description = [sample,gene,var,var,description["product"]]
                    flattened_vars.append(variant_description)
    return flattened_vars

def get_gb_coords(gb_file, gp_list):
    '''
    DEPRECATED
    Function to take a main genbank 'coordinates' file and a list of any genpept 
    files to replace main gene coordinates with. E.g. in Sars-CoV-2 ORF1ab gene 
    produces multiple mature peptides and may want variants to be annotated in 
    reference to position within mature peptide rather than overall gene.  
    '''
    genbank = SeqIO.read(open(gb_file,"r"), "genbank")
    genpepts = {}
    for item in gp_list:
        genpept = SeqIO.read(open(item, "r"), "genbank")
        for feature in genpept.features:
            if feature.type == "Protein":
                genpepts[feature.qualifiers["product"][0]] = genpept
    coords = {}
    for feature in genbank.features:
        if feature.type == "CDS": #the coding sequences define the gene products
            if feature.qualifiers["product"][0] in genpepts.keys():
                if feature.qualifiers["gene"][0] in coords:
                    for mat_pep in genpepts[feature.qualifiers["product"][0]].features:
                        if mat_pep.type == "mat_peptide":
                            coords[feature.qualifiers["gene"][0]]["multi-product"][mat_pep.qualifiers["product"][0]] = {"sub-position": [int(mat_pep.location.start) ,int(mat_pep.location.end -1)]}
                            coords[feature.qualifiers["gene"][0]]["product"] = feature.qualifiers["gene"][0]
                else:
                    coords[feature.qualifiers["gene"][0]]= {}
                    coords[feature.qualifiers["gene"][0]]["multi-product"] = {}
                    coords[feature.qualifiers["gene"][0]]["product"] = {}
                    for mat_pep in genpepts[feature.qualifiers["product"][0]].features:
                        if mat_pep.type == "mat_peptide":
                            coords[feature.qualifiers["gene"][0]]["multi-product"][mat_pep.qualifiers["product"][0]] = {"sub-position": [int(mat_pep.location.start) ,int(mat_pep.location.end -1)]}
                            coords[feature.qualifiers["gene"][0]]["product"] = feature.qualifiers["gene"][0]
            else:
                coords[feature.qualifiers["gene"][0]] = {"product" : feature.qualifiers["product"][0]}

    return coords
    

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--genpept', metavar='', nargs='+',
        help='secondary genpept files for matured peptide annotation e.g. --genpept data/*.gp')
    parser.add_argument('--json', metavar='output.variants.json', type=str,
        help='Destination file for corrected variant list (JSON)')
    parser.add_argument('--csv', metavar = 'sample_summary.csv', type = str)
    parser.add_argument('input_csv', metavar='input.variants.csv', type=str,
        help='Input variant csv file')
    parser.add_argument('genbank', metavar='reference.gb', type=str,
        help='primary genbank file for annotation')

    args = parser.parse_args()
    start_msg = f'''
    Running with the following parameters:
                 Input file: {args.input_csv}
                Output file: {args.json}
        Output CSV Location: {args.csv}
               Genbank file: {args.genbank}
            Genpept file(s): {args.genpept}
    '''
    print(start_msg)
    
    coords = get_gb_coords(args.genbank, args.genpept)
    samples = variants_parser(args.input_csv)
    #passing a deepcopy of the nested dictionary structure to prevent function modifying the mutable input dictionary. Not sure when this would be a problem but better to keep them separate
    corrected_samples = get_variant_info(copy.deepcopy(samples),copy.deepcopy(coords)) 
    if args.json:
        with open(args.json, 'w') as outfile:
            json.dump(corrected_samples, outfile)
    if args.csv:
        flattened_variants = flatten_variants(corrected_samples)
        with open(args.csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(flattened_variants)

if __name__ == "__main__":
    main()