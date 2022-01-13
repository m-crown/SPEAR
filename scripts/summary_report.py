#!/usr/bin/env python

from pandas.io.formats.format import SeriesFormatter
from Bio.SeqUtils import seq1
from Bio import SeqIO
import pandas as pd
import argparse
from pathlib import Path
import numpy as np
from summarise_snpeff import parse_vcf, write_vcf
import csv
import re

def main():