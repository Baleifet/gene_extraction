'''
a python script to extract protein sequences from genbank files given a search string for a specific annotation
    arguments:
        argv[1]:
        argv[2]:
        argv[3]:
        argv[4]:
   
    functions:
        x
        x
        x
'''
from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
 
try: 
    feature_type = argv[1]
    search_string = argv[2]
    ext = argv[3] ## implement this!
    directory = argv[4]
except:
    print("...")
    exit()

def parse_genbank_files(file_name,directory):
    #looks through all files in a directory and makes a temporary lsit of records
    records = []
    file_path = os.path.join(directory, file_name)
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            records.append(record)
    return records
 
def proteinsquarry(records, annotation=None,featuretype=None):
    #searches thorough all the records and returns those annotated as wanted of a specific feature_type
    translations={}
    annotation = annotation.lower() #if annotation is not None and isinstance(annotation, str) else annotation
 
    if featuretype is not None and annotation is not None:
        for i in range(len(records)):
                translations.update({f.qualifiers.get("locus_tag",[None])[0]:[
                    f.qualifiers.get("product",[None])[0],
                    f.qualifiers.get("translation",[None])[0]]
                    for f in records[i].features if f.type==featuretype and f.qualifiers.get("product") is not None and annotation in f.qualifiers.get("product")[0].lower()})
                #for f in records[i].features:
                #    print(annotation)
                #    if f.type==featuretype and f.qualifiers.get("product") is not None and annotation in f.qualifiers.get("product")[0].lower():
                #        print('wow')
 
    elif featuretype is not None and annotation==None:
        #add annotations
        for i in range(len(records)):
                translations.update({f.qualifiers.get("locus_tag",[None])[0]:[
                    f.qualifiers.get("product",[None])[0],
                    f.qualifiers.get("translation",[None])[0]]
                    for f in records[i].features if f.type==featuretype})
 
    elif annotation is not None and featuretype==None:
        #add featuretypes
        for i in range(len(records)):
                translations.update({f.qualifiers.get("locus_tag",[None])[0]:[
                    f.qualifiers.get("product",[None])[0],
                    f.qualifiers.get("translation",[None])[0]]
                    for f in records[i].features if (f.qualifiers.get("product") is not None and annotation in f.qualifiers.get("product")[0].lower())})
    return translations
 
def FASTAexport(sequences, file_name):
    #exports all sequences in a FASTA format
    #'sequences' must be a dictinoary
    to_export = []
    for name, sequence in sequences.items():
        seq_record = SeqRecord(Seq(sequence[1]), id=name, description=" - %s"%sequence[0])
        to_export.append(seq_record)
    SeqIO.write(to_export, file_name, "fasta")
 
#Main Routine
translations_combined = {}

for file_name in os.listdir(directory):
    try:
        if file_name.endswith(".gbk"): #!!
            records = parse_genbank_files(file_name,directory)
            print('found', len(records), search_string, 'in', file_name)
            translations_combined.update(proteinsquarry(records, search_string, feature_type))
    except Exception as error:
        print("There was an error reading:", file_name, str(error))
FASTAexport(translations_combined, f'{os.path.basename(directory)}_{feature_type}_{search_string}.faa')