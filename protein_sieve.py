'''
a python script to extract protein sequences from genbank files given necessarily a directory folder to parse from, and various arguments such as: strings to search for in the name, strings to avoid, a feature type to select for, feature types to avoid, minimum length of the translation.
usage example: --dir DIR --search integrase, recombinase, dehydrogenase --negsearch tyrosine --featype CDS --length 100
    will return all CDS features that contain the search words in the annotation name but not the word tyrosine, from all .gbk files in the directory DIR and whose translation is at least 100 amino acids long.
usage example: --search integrase, --dir .  --negsearch PFAM, PFAM_domain
    will return all integrases in all .gbk files in the directory DIR that were NOT annotated as 'PFAM' or 'PFAM_domain'.

    arguments:
        --featype     : feature type (e.g. CDS, PFAM_fomain, etc)
        --search      : search one or more strings (e.g. integrase, recombinase, dehydrogenase)
        --search_all  : search all strings present at the same time in the name, non-ordered 
        --ext         : files extension (e.g. .gbk)
        --dir         : directory to parse from
        --length      : minimum length of the feature's translation
        --negsearch   : strings to avoid
        --negfeatype  : feature type to avoid
   
    functions:
        make_list(arguments)
        parse_genbank_files(file_name,directory)
        add_protein(f,qualifier,translations,file_name)
        get_qualifier(f)
        proteinsquarry(records, file_name, search=None, negsearch=None, feature_type=None, negfeature_type=None)
        export_sequences(sequences,file_name,map)
'''

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time

parser = argparse.ArgumentParser()
parser.add_argument("--dir", help="Directory containing files to parse through", required=True)
parser.add_argument("--search", nargs='*', help="A string or any among a group of strings to search for", default=None)
parser.add_argument("--search_all", nargs='*', help="Search all strings among a list of strings", default=None)
parser.add_argument("--featype", nargs='*', help="Feature type to search for", default=None) 
parser.add_argument("--ext", help="File extension to filter results", default='.gbk')
parser.add_argument("--length", help="Minimum length of features to consider", type=int, default=None)
parser.add_argument("--negsearch", nargs='*', help="A string or a list of strings to exclude from the research", default=None) #pseudogene etc
parser.add_argument("--negfeatype", nargs='*', help="Feature type to exclude from the research", default=None)


def make_list(arguments):
    if '[' in arguments or '(' in arguments:
        arguments = arguments.replace('[','').replace(']','').replace('(','').replace(')','')
    if ',' in arguments:
        arguments = arguments.split(',')
        arguments = [argument.replace("'","").replace('"','').replace(' ','') for argument in arguments]
        return [argument for argument in arguments if argument != '']
    else:
        arguments=[arguments.replace("'","").replace('"','').replace(' ','')]
        return arguments

args = parser.parse_args()
directory = str(args.dir)
feature_type = make_list(str(args.featype))             if args.featype is not None else None
search = make_list(str(args.search))                    if args.search is not None else None
search_all = make_list(str(args.search_all))            if args.search_all is not None else None
ext = str(args.ext)                                     if args.ext is not None else None
min_length = int(args.length)                           if args.length is not None else None
negsearch = make_list(str(args.negsearch))              if args.negsearch is not None else None
negfeature_type = make_list(str(args.negfeatype))       if args.negfeatype is not None else None

start_time = time.time()

def parse_genbank_files(file_name,directory):
    # looks through all files in a directory and makes a temporary lsit of records
    records = []
    file_path = os.path.join(directory, file_name)
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "genbank"): #???
            records.append(record)
    return records

def add_protein(f,qualifier,translations,file_name):
    # adds the protein (f) to the 'translations' dictionary
    try:
        if min_length is not None and len(f.qualifiers.get("translation",[None])[0]) > min_length:
            translations.update({f.qualifiers.get("locus_tag",[None])[0]:
                                    [f.qualifiers.get(qualifier,[None])[0],
                                    f.qualifiers.get("translation",[None])[0],
                                    len(f.qualifiers.get("translation",[None])[0]),
                                    f.type,
                                    f.location,
                                    file_name]})
        elif min_length is None:
            translations.update({f.qualifiers.get("locus_tag",[None])[0]:
                                    [f.qualifiers.get(qualifier,[None])[0],
                                    f.qualifiers.get("translation",[None])[0],
                                    len(f.qualifiers.get("translation",[None])[0]),
                                    f.type,
                                    f.location,
                                    file_name]})
        return translations
    except TypeError:
        pass

def get_qualifier(f):
    # this function defines which name to give to each feature_type depending on which of the dictionary values appears in the feature.qualifiers. checked in a priority order, can be changed as needed
    qualifier_names = {
        'CDS' : ['product'],
        'CDS_motif' : ['domain_id', 'aSTool'],
        'PFAM' : ['description','domain_id'],
        'PFAM_domain' : ['description','domain_id'],
        'aSDomain' : ['description','domain_id']
    }
    qualifier = None
    for name in qualifier_names.get(str(f.type),[None]):
        if f.qualifiers.get(name,[None])[0] is not None:
            qualifier = name
    return qualifier

def proteinsquarry(records, file_name, search=None, search_all=None, negsearch=None, feature_type=None, negfeature_type=None):
    translations={}
    # for every protein
    for i in range(len(records)):
        for f in records[i].features:
            # the feature can either be of the desired type or not be of the undesired type
            if (feature_type is not None and any(feature_type_string == str(f.type) for feature_type_string in feature_type)) or (negfeature_type is not None and not any(negfeature_type_string == str(f.type) for negfeature_type_string in negfeature_type)):
                if "translation" in f.qualifiers:
                    for qualifier in f.qualifiers:
                        # search through all qualifiers the ones that contain ANY of the strings contained in --search arg, but not in the --negsearch if this is specified, and only if the feature has a protein sequence ("translation")
                        if search is not None and any(search_string.lower() in f.qualifiers.get(qualifier)[0].lower() for search_string in search):
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name)
                        # search through all qualifiers the ones that contain ALL of the strings contained in --search arg, but not in the --negsearch if this is specified, and only if the feature has a protein sequence ("translation")
                        if search_all is not None and all(search_all_string.lower() in f.qualifiers.get(qualifier)[0].lower().split(' ') for search_all_string in search_all):
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name)
                        # if both search and search_all are not specified, it will output all the features that are of the desired type but don't contain any of the --negsearch strings
                        if search is None and search_all is None:
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name)
            elif feature_type is None and negfeature_type is None:
                # if there's neither a desired nor an undesired feature_type, loops through all qualifiers searching for search, just as above, but without feature_type constrain
                if "translation" in f.qualifiers:
                    for qualifier in f.qualifiers:
                        if search is not None and any(search_string.lower() in f.qualifiers.get(qualifier)[0].lower() for search_string in search):
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name)
                        if search_all is not None and all(search_all_string.lower() in f.qualifiers.get(qualifier)[0].lower().split(' ') for search_all_string in search_all):
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name)
                        if search is None and search_all is None:
                            if negsearch is not None and not any(negsearch_string.lower() in f.qualifiers.get(qualifier)[0].lower() for negsearch_string in negsearch):
                                add_protein(f,get_qualifier(f),translations,file_name)
                            elif negsearch is None:
                                add_protein(f,get_qualifier(f),translations,file_name) # None of the arguments have been provided, this call will add all the features from all the files to the output
    return translations

def export_sequences(sequences,title,map):
    # exports an output summary in two .csv files and all sequences in a FASTA file
    os.makedirs(title, exist_ok=True)

    # setup the two dataframe output
    df1 = pd.DataFrame(list(map.items()),columns=['file_name','number of proteins'])
    df2 = pd.DataFrame(sequences)
    df1.to_csv(os.path.join(title, f'{title}_map.csv'), index=True)
    first_key, first_value = next(iter(df2.items()))
    df2=df2.transpose()
    df2.columns = ['qualifier name','translation','length (translation)','feature type','location','source file']
    df2.to_csv(os.path.join(title, f'{title}_details.csv'), index=True)

    to_export=[]
    if type(sequences) == dict:
        for name, sequence in sequences.items():
            seq_record = SeqRecord(Seq(sequence[1]), id=name, description='')
            to_export.append(seq_record)
        SeqIO.write(to_export, os.path.join(title, f'{title}_sequences.faa'), 'fasta')
    else:
        print("'sequences_dict' must be a dictinoary")
 
 
if directory == '.':
    real_directory = os.path.basename(os.getcwd())
else:
    real_directory = os.path.basename(directory)

def generate_title(args):
    title_parts = [real_directory]
    if args.search is not None:
        title_parts.append(f'search({str(search)[1:-1]})')
    if args.search_all is not None:
        title_parts.append(f'search_all({str(search_all)[1:-1]})')
    if args.featype is not None:
        title_parts.append(f'featype({str(feature_type)[1:-1]})')
    if args.length is not None:
        title_parts.append(f'length({args.length})')
    if args.negsearch is not None:
        title_parts.append(f'negsearch({str(negsearch)[1:-1]})')
    if args.negfeatype is not None:
        title_parts.append(f'negfeatype({str(negfeature_type)[1:-1]})')
    title = '_'.join(title_parts)
    return title
title = generate_title(args)
 
# Main Routine
translations_combined = {}
map = {}

for file_name in os.listdir(directory):
    try:
        if file_name.endswith(ext): 
            records = parse_genbank_files(file_name,directory)
            translations = proteinsquarry(records, file_name, search, search_all, negsearch, feature_type, negfeature_type)
            translations_combined.update(translations)
            map.update({file_name:str(len(translations))})
            print('found', len(translations),'results in', file_name)
        else:
            continue
    except Exception as error:
        print("There was an error reading:", file_name, str(error))
print("\nfound", len(translations_combined), "total results in", real_directory)
print(f"\nsaved as {title} in {os.getcwd()} \n")

export_sequences(translations_combined, generate_title(args), map)

end_time = time.time()
total_wall_time = end_time - start_time
print(f"Total Wall Time: {total_wall_time} seconds\n")