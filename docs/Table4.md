## SPEAR Score Summary
The table below outlines the fields in `spear_score_summary.tsv` which is produced per analysis run and has one row for each sample.

| Column ID | Description |
| --------- | ----------- |
| `sample_id` | Input sample ID, taken from input file header in `.fa`, `.aln` and sample col in `.vcf` |
| `total_variants` | Total number of genomic variants at nucleotide level  |
| `total_residue_variants` | The total number of amino acid changes, here, note insertions are counted as a single event, each deleted residue is counted |
| `consequence_type_variants` | Comma separated list of `consequence:count` present in sample, summary of values in `consequence_type` from Table 2 "
| `region_residues` | Amino acid changes as comma separated list summarised per `product:region:count` |
| `domain_residues` | Amino acid changes as comma separated list summarised per `product:domain:count`, domain defitions can be found in [`/docs`](docs/) |
| `ACE2_contact_counts` | Total number of mutated amino acids involved in ACE2 contacts |
| `ACE2_contact_score` | Sum of contact scores, salt bridges:3, hydrogen-bonds:2, generic residue contact:1 |
| `trimer_contact_counts` | Total number of mutated amino acids involved in trimer interface contacts, this is specifically for Spike |
| `trimer_contact_score` | Sum of contact scores for Spike trimer interface, slat bridges:3, hydrogen-bonds:2, generic residue contact:1 |
| `barns_class_variants` | Total number of mutated residues per barns class epitope summarised as comma separated list, as `class:count` additional residues not formally part of epitope but with strong Deep Mutational Scanning (DMS) evidence as having an impact on binding of mAbs of that class are flagged with `*`. |
| `bloom_ACE2_sum` | Sum of ACE2 binding value Δlog10(KD,app) relative to the "wild-type" (WT) across all mutated residues, data from Starr _et al_.[3] Higher positive values mean binding is stronger than WT, negative values mean binding is weaker then WT |
| `bloom_ACE2_max` | Highest scoring residue for ACE2 value out of all mutated residues |
| `bloom_ACE2_min` | Lowest scoring residue for ACE2 value out of all mutated residues |
| `VDS_sum` | Sum of all Vibrational Difference Scores (VDS), positive VDS values suggests mutation stabilises the open state of Spike and/or makes the closed state more flexible, favouring the open conformation relative to the WT. Negative values suggest mutation favours the closed state more than WT. Data from Teruel _et al_.[4] |
| `VDS_max`| Highest scoring residue for VDS value out of all mutated residues |
| `VDS_lowest`| Lowest scoring residue for VDS value out of all mutated residues |
| `serum_escape_sum` | Sum of mean residue specific serum escape score from 7 individuals in Greaney _et al_.[5] across all mutations, larger values between indicate more escape |
| `serum_escape_max` | Highest scoring residue for serum escape out of all mutated residues |
| `serum_escape_min`	| Lowest scoring residue for serum escape out of all mutated residues| 
| `mAb_escape_all_classes_sum` | Sum of mean residue specific mAb escape scores for all mutated residues in sample from 26 mAbs data taken from Dong _et al_.[6], SARS-CoV-2-RBD\_MAP\_COV2-2955[7], Greaney _et al_.[8], Starr _et al_.[9], Starr _et al_.[10], Starr _et al_.[11] Tortorici _et al_.[12] |
| `mAb_escape_all_classes_max` | Highest scoring residue for mAb escape out of all mutated residues |
| `mAb_escape_all_classes_min` | Lowest scoring residue for mAb escape out of all mutated residues |
| `mAb_escape_class_1_sum` | As above, sum of mean residue specific mAb escape score for all mutated residues from class 1 mAbs only, only applied to residues in Barns class 1 epitope |
| `mAb_escape_class_1_max` | Highest scoring residue for class 1 mAb escape out of all mutated residues |
| `mAb_escape_class_1_min` | Lowest scoring residue for class 1 mAb escape out of all mutated residues |
| `mAb_escape_class_2_sum` | As above, sum of mean residue specific mAb escape score for all mutated residues from class 2 mAbs only, only applied to residues in Barns class 2 epitope |
| `mAb_escape_class_2_max` | Highest scoring residue for class 2 mAb escape out of all mutated residues |
| `mAb_escape_class_2_min` | Lowest scoring residue for class 1 mAb escape out of all mutated residues |
| `mAb_escape_class_3_sum` | As above, sum of mean residue specific mAb escape score for all mutated residues from class 3 mAbs only, only applied to residues in Barns class 3 epitope | 
| `mAb_escape_class_3_max` | Highest scoring residue for class 3 mAb escape out of all mutated residues |
| `mAb_escape_class_3_min` | Lowest scoring residue for class 1 mAb escape out of all mutated residues | 
| `mAb_escape_class_4_sum` | As above, sum of mean residue specific mAb escape score for all mutated residues from class 4 mAbs only, only applied to residues in Barns class 4 epitope | 
| `mAb_escape_class_4_max` | Highest scoring residue for class 4 mAb escape out of all mutated residues | 
| `mAb_escape_class_4_min` | Lowest scoring residue for class 4 mAb escape out of all mutated residues | 
| `BEC_RES_sum` | Sum of all per residue Bloom Escape Calculator Residue Escape scores |
| `BEC_RES_max` | Highest BEC scoring residue |
| `BEC_RES_min` | Lowest BEC scoring residue | 
| `BEC_EF_sample` | Bloom Escape Calculator Escape factor as `BEC_EF` but calculated using the full compliment of changed residues in the whole sample |

**Table 4**. `spear_score_summary.tsv` fields.