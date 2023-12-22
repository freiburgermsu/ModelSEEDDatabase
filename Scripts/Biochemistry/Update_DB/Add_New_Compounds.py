#!/usr/bin/env python
import os, sys, re, copy
import argparse
from csv import DictReader
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument('compounds_file', help="Compounds File")
parser.add_argument('database', help="Biochemistry database of origin")
parser.add_argument("-s", dest='save_file', action='store_true')
parser.add_argument("-r", dest='report_file', action='store_true')
parser.add_argument("-no", dest='names_only', action='store_true')
args = parser.parse_args()

if(os.path.isfile(args.compounds_file) is False):
    print("Cannot find file: "+args.compounds_file)
    sys.exit()

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

original_name_dict = compounds_helper.loadNames()
searchnames_dict = dict()
all_names_dict = dict()
new_name_count = dict()
for msid in sorted(original_name_dict):
    for name in original_name_dict[msid]:
        all_names_dict[name]=1

        searchnames = compounds_helper.searchname(name)
        for searchname in searchnames:
            #Avoid redundancy where possible
            if(searchname not in searchnames_dict):
                searchnames_dict[searchname]=msid

original_alias_dict=compounds_helper.loadMSAliases()
new_alias_count = dict()
source_alias_dict = dict()
for msid in original_alias_dict:
    for source in original_alias_dict[msid]:
        if(source not in source_alias_dict):
            source_alias_dict[source]=dict()

        for alias in original_alias_dict[msid][source]:
            if(alias not in source_alias_dict[source]):
                source_alias_dict[source][alias]=list()
            source_alias_dict[source][alias].append(msid)
            
unique_structures_dict = compounds_helper.loadStructures(["InChI","InChIKey","SMILE"],["ModelSEED"])
all_structures_dict = compounds_helper.loadStructures(["InChI","InChIKey","SMILE"],["ModelSEED"],False)
compiled_structures_dict = dict()
for msid in unique_structures_dict:
    for structure_format in ["InChI","InChIKey","SMILE"]:
        if(structure_format in unique_structures_dict[msid]):
            for structure in unique_structures_dict[msid][structure_format]:
                struct_list = [structure]
                if(structure_format == 'InChIKey'):
                    struct_list.append('-'.join(structure.split('-')[0:2]))
                    #struct_list.append('-'.join(structure.split('-')[0:1]))

                for struct in struct_list:
                    if(struct not in compiled_structures_dict):
                        compiled_structures_dict[struct]={'Unique':[],'Charged':[],'Original':[]}
                    compiled_structures_dict[struct]['Unique'].append(msid)

# We keep a list of structures for which we haven't curated the conflicts
# we can still attempt to match these, to give us an idea
for msid in all_structures_dict:
    for structure_format in ["InChI","InChIKey","SMILE"]:
        if(structure_format in all_structures_dict[msid]):
            for structure in all_structures_dict[msid][structure_format]:
                struct_list = [structure]
                struct_stage = all_structures_dict[msid][structure_format][structure]['type']

                if(structure_format == 'InChIKey'):
                    struct_list.append('-'.join(structure.split('-')[0:2]))
                    #struct_list.append('-'.join(structure.split('-')[0:1]))

                for struct in struct_list:
                    if(struct not in compiled_structures_dict):
                        compiled_structures_dict[struct]={'Unique':[],'Charged':[],'Original':[]}
                    compiled_structures_dict[struct][struct_stage].append(msid)

#Find last identifier and increment
last_identifier = list(sorted(compounds_dict))[-1]
identifier_count = int(re.sub('^cpd','',last_identifier))

Default_Cpd = OrderedDict({ "id":"cpd00000","name":"null","abbreviation":"null","aliases":"null",
                             "formula":"null","mass":10000000,"charge":0,
                             "deltag":10000000.0,"deltagerr":10000000.0,"pka":"","pkb":"",
                             "inchikey":"","smiles":"",
                             "is_cofactor":0,"is_core":0,"is_obsolete":0,
                             "abstract_compound":"null","comprised_of":"null","linked_compound":"null",
                             "notes":[],"source":"" })
New_Cpd_Count=dict()
matched_cpds_dict=dict()
Headers=list()
with open(args.compounds_file) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        cpd=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            cpd[Headers[i].lower()]=array[i]

        matched_cpd = {'msid':None,'string':None,'format':None,'stage':None}
        cpd_has_structure=False

        #First check that the Alias doesn't already exist
        if(args.database in source_alias_dict and cpd['id'] in source_alias_dict[args.database]):
            msids = list(sorted(source_alias_dict[args.database][cpd['id']]))
            if(len(msids)>0):
                matched_cpd['msid']=msids[0]
                matched_cpd['format']='id'
                matched_cpd['string']=cpd['id']
                matched_cpd['stage']='na'

        #Then check that the Structure doesn't already exist, first as InChI, then as InChiKey, then as SMILES
        if(matched_cpd['msid'] is None and args.names_only is False):

            for struct_format in ['inchi','inchikey','smile','smiles']:
                if(struct_format in cpd):
                    cpd_has_structure=True
                    structure_list = [cpd[struct_format]]
                    if(struct_format == 'inchikey'):
                        structure_list.append('-'.join(cpd[struct_format].split('-')[0:2]))
                        structure_list.append('-'.join(cpd[struct_format].split('-')[0:1]))

                    for structure in structure_list:
                        
                        if(matched_cpd['msid'] is not None):
                            break

                        if(structure in compiled_structures_dict):

                            for struct_stage in ['Unique','Charged','Original']:
                                
                                if(matched_cpd['msid'] is not None):
                                    break

                                msids = list(sorted(compiled_structures_dict[structure][struct_stage]))
                                if(len(msids)>0):
                                    matched_cpd['msid']=msids[0]
                                    matched_cpd['format']=struct_format
                                    matched_cpd['string']=structure
                                    matched_cpd['stage']=struct_stage

        #Then check that the Name doesn't already exist
        if(matched_cpd['msid'] is None):
            msids=dict()
            for name in cpd['names'].split('|'):
                searchnames = compounds_helper.searchname(name)
                for searchname in searchnames:
                    if(searchname in searchnames_dict):
                        # If we're searching names, that means the structures didn't match
                        # But, if the new cpd has a structure, it can't match any compound
                        # that already has a structure, hence this condition
                        if(cpd_has_structure is True and \
                            searchnames_dict[searchname] in all_structures_dict):  
                            continue

                        if(searchnames_dict[searchname] not in msids):
                            msids[searchnames_dict[searchname]]=name
                            break
                        
            msids_list=list(sorted(msids))
            if(len(msids_list)>0):
                matched_cpd['msid']=msids_list[0]
                matched_cpd['format']='name'
                matched_cpd['string']=msids[msids_list[0]]
                matched_cpd['stage']='na'

        if(cpd['id'] not in matched_cpds_dict):
            matched_cpds_dict[cpd['id']]=list()
        matched_cpds_dict[cpd['id']].append(matched_cpd)

        if(matched_cpd['msid'] is not None):

            #Regardless of match-type, add new names
            #NB at this point, names shouldn't match _anything_ already in the database
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['names'].split('|'):
                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(matched_cpd['msid'] not in original_name_dict):
                        original_name_dict[matched_cpd['msid']]=list()
                    original_name_dict[matched_cpd['msid']].append(name)
                    all_names_dict[name]=1
                    new_name_count[matched_cpd['msid']]=1
                
            #if matching structure or name, add ID to aliases
            if(matched_cpd['format'] != 'id'):
                if(matched_cpd['msid'] not in original_alias_dict):
                    original_alias_dict[matched_cpd['msid']]=dict()
                if(matched_cpd['msid'] in original_alias_dict and args.database not in original_alias_dict[matched_cpd['msid']]):
                    original_alias_dict[matched_cpd['msid']][args.database]=list()
                original_alias_dict[matched_cpd['msid']][args.database].append(cpd['id'])
                new_alias_count[matched_cpd['msid']]=1

            #Update source type
            compounds_dict[matched_cpd['msid']]['source']='Primary Database'

        elif(args.save_file is True):

            #New Compound!
            #Generate new identifier
            identifier_count+=1
            new_identifier = 'cpd'+str(identifier_count)

            new_cpd = copy.deepcopy(Default_Cpd)
            new_cpd['id']=new_identifier
            if('mass' in cpd):
                new_cpd['mass']=float(cpd['mass'])
            if('charge' in cpd):
                new_cpd['charge']=int(cpd['charge'])
            if('formula' in cpd):
                new_cpd['formula']=cpd['formula']

            #Add new identifier with original ID as alias
            original_alias_dict[new_cpd['id']]={args.database:[cpd['id']]}
            new_alias_count[new_cpd['id']]=1

            #Add new names
            #Names are saved separately as part of the aliases at the end of the script
            for name in cpd['names'].split('|'):
                if(new_cpd['name']=='null'):
                    new_cpd['name']=name
                    new_cpd['abbreviation']=name

                if(name not in all_names_dict):
                    #Possible for there to be no names in biochemistry?
                    if(new_cpd['id'] not in original_name_dict):
                        original_name_dict[new_cpd['id']]=list()
                    original_name_dict[new_cpd['id']].append(name)
                    all_names_dict[name]=1
                    new_name_count[new_cpd['id']]=1

            #If no names at all
            if(new_cpd['name']=='null'):
                new_cpd['name']=cpd['id']
                new_cpd['abbreviation']=cpd['id']

            #Add source type
            new_cpd['source']='Primary Database'
            compounds_dict[new_cpd['id']]=new_cpd
            New_Cpd_Count[new_cpd['id']]=1

#Here, for matches, re-write names and aliases
print("Compounds matched via:")
matched_src_dict=dict()
for oc in sorted(matched_cpds_dict):
    for matched_cpd in matched_cpds_dict[oc]:
        if(matched_cpd['format'] not in matched_src_dict):
            matched_src_dict[matched_cpd['format']]=list()
        matched_src_dict[matched_cpd['format']].append(matched_cpd['msid'])

for src in matched_src_dict:
    if(src is not None):
        match_dict = dict()
        for match in matched_src_dict[src]:
            match_dict[match]=1
        print("\t"+src+" "+str(len(matched_src_dict[src]))+" matched to "+str(len(match_dict.keys())))
if(None in matched_src_dict):
    print("\t"+str(len(matched_src_dict[None]))+" not matched to any ModelSEED compounds")
        
if(args.report_file is True):
    file_stub = '.'.join(args.compounds_file.split('.')[0:-1])
    report_file = file_stub+'.rpt'
    print("Saving report to file: "+report_file)
    with open(report_file,'w') as rfh:
        for oc in sorted(matched_cpds_dict):
            for matched_cpd in matched_cpds_dict[oc]:
                rfh.write(oc+'\t'+str(matched_cpd['msid'])+'\t'+str(matched_cpd['format']))
                rfh.write('\t'+str(matched_cpd['string'])+'\t'+str(matched_cpd['stage'])+'\n')
        
if(args.save_file is True):
    print("Saving additional names for "+str(len(new_name_count))+" compounds")
    compounds_helper.saveNames(original_name_dict)
    print("Saving additional "+args.database+" aliases for "+str(len(new_alias_count))+" compounds")
    compounds_helper.saveAliases(original_alias_dict)
    print("Saving "+str(len(New_Cpd_Count))+" new compounds from "+args.database)
    compounds_helper.saveCompounds(compounds_dict)

#Scripts to run afterwards
#./Merge_Formulas.py
#./Update_Compound_Aliases.py
#../Structures/List_ModelSEED_Structures.py
#../Structures/Update_Compound_Structures_Formulas_Charge.py
#./Rebalance_Reactions.py