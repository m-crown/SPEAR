#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Matt Crown 2021"""


from Bio import SeqIO
import csv
import copy
import argparse
import json

def variants_parser(variants_list):
    '''
    Takes the ouput CSV from gofasta sam variants and stores samples in a nested dictionary containing lists of mutations for each gene.
    '''
    samples = {}
    with open(variants_list, 'r') as file:
        csvreader = csv.reader(file)
        header = next(csvreader)
        for sample in csvreader:
            samples[sample[0]] = {}
            for var in sample[1].split("|"): #split out pipe delimited variants
                gene, pos = var.split(":") #split variant call into gene and position
                if gene in samples[sample[0]].keys(): #check if the gene already exists in nested dictionary keys
                    if pos not in samples[sample[0]][gene].keys(): #deduplicate the output from gofasta sam variants - particularly needed for orf1ab variants which are in orf1a and orf1ab
                        samples[sample[0]][gene][pos] = {} #append position to exiting list for this gene
                else:
                    samples[sample[0]][gene] = {pos: {}} #if gene doesnt already exist in nested dictionary create the key and set values to be a list including the position
    return samples

def get_variant_info(samples, coords):
    '''
    Function takes a list of variants with polyprotein positions and converts them to matured peptide positions; also returns name of the matured peptide the variant is in
    '''
    for sample in samples.values():        
        for gene_name, variants_list in sample.items():
            if gene_name != "synSNP":
                for variant_name, variant_info in variants_list.items():
                    variant_info["product"] = coords[gene_name]["product"]
                    if "multi-product" in coords[gene_name]:
                        for k,v in coords[gene_name]["multi-product"].items():
                            if v["sub-position"][0] <= int(variant_name[1:-1]) -1 <= v["sub-position"][1]:
                                variant_info["multi-product"] = k
                                variant_info["sub-position"] = f'{variant_name[0]}{int(variant_name[1:-1]) - v["sub-position"][0]}{variant_name[-1]}'
                                break
    return samples

#NEED TO WRITE A FUNCTION TO PARSE THE OUTPUT AND RETRIEVE VARIANTS IN CSV FORMAT? OR JUST HAVE IT OUTPUT A CSV VERSION TOO 

def main():
    parser = argparse.ArgumentParser(description='Encrypts user defined fields in JSON file with supplied encryption key')
    parser.add_argument('-i', '--input_csv', metavar='', type=str,
                    help='Input variant csv file')
    parser.add_argument('-o', '--output_file', metavar='', type=str,
                    help='Destination file for updated variant list')
    parser.add_argument('-d', '--data', metavar='', type=str,
                    help='location of genbank/genpept files')

    args = parser.parse_args()

    print("Running with the following parameters:\nInput file: ", args.input_csv, "\nOutput file: ", args.output_file, "\n", "Data Loc: ", args.data, "\n")

    # #check what the input file type is (csv or json) - this is probably not necessary for dehashing a report but keeping in just in case its needed
    # if args.input_file.endswith('.json'):
    #     input_file_type = "JSON"
    #     with open(args.input_file) as json_file:
    #         data = json.load(json_file)



    orf1ab_gb = SeqIO.read(open(f'{args.data}/ORF1ab.gp',"r"), "genbank")
    sarscov2_gb = SeqIO.read(open(f'{args.data}/NC_045512.2.gb', "r"), "genbank")
    sarscov2_coords = {}
    for feature in sarscov2_gb.features:
        if feature.type == "CDS":
            if feature.qualifiers["gene"][0] == "ORF1ab":
                sarscov2_coords["ORF1ab"]= {}
                sarscov2_coords["ORF1ab"]["multi-product"] = {}
                sarscov2_coords["ORF1ab"]["product"] = {}
                for mat_pep in orf1ab_gb.features:
                    if mat_pep.type == "mat_peptide":
                        sarscov2_coords["ORF1ab"]["multi-product"][mat_pep.qualifiers["product"][0]] = {"sub-position": [int(mat_pep.location.start) ,int(mat_pep.location.end -1)]}
                        sarscov2_coords["ORF1ab"]["product"] = feature.qualifiers["gene"][0]
            else:
                sarscov2_coords[feature.qualifiers["gene"][0]] = {"product" : feature.qualifiers["product"][0]}

    go_fasta_samples = variants_parser("variants.csv")
    samples = get_variant_info(copy.deepcopy(go_fasta_samples),copy.deepcopy(sarscov2_coords)) #passing a deepcopy of the nested dictionary structure to prevent function modifying the mutable input dictionary. Not sure when this would be a problem but better to keep them separate
    with open(args.output_file, 'w') as outfile:
        json.dump(samples, outfile)

if __name__ == "__main__":
    main()