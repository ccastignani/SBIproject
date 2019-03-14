import os
import sys
import tarfile
from Bio.PDB import *
import re

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def superimposition(struc_ref,struc_sample,chain_ref,chain_sample,current_chain):
	parser = PDBParser()
	if type(struc_ref) is str:
		ref_structure = parser.get_structure("reference", struc_ref)
	else: 
		ref_structure = struc_ref
	sample_structure = parser.get_structure("sample", struc_sample)

	ref_model    = ref_structure[0]
	sample_model = sample_structure[0]

	ref_atoms = []
	sample_atoms = []

	for ref_chain in ref_model:
		if ref_chain.get_id() == chain_ref:
			#Iterate of all residues in each model in order to find proper atoms
			for ref_res in ref_chain:
				ref_atoms.append(ref_res['CA'])

			# Do the same for the sample structure
			for sample_chain in sample_model:
				if sample_chain.get_id() == chain_sample:
					for sample_res in sample_chain:
						sample_atoms.append(sample_res['CA'])

	super_imposer = Superimposer()
	super_imposer.set_atoms(ref_atoms, sample_atoms)
	super_imposer.apply(sample_model.get_atoms())

	rotmatrix = super_imposer.rotran[0]
	translation = super_imposer.rotran[1]

	#Create new chain for the sueprimposition result
	new_chain=Chain.Chain(chr(current_chain))
	
	ref_structure[0].add(new_chain)
	for chain in sample_model:
		if chain.get_id()!=chain_sample:
			for residue in chain.get_residues():
				res=residue.get_id()[1]
				ref_structure[0][chr(current_chain)].add(residue)
				for atom in ref_structure[0][chr(current_chain)][res].get_atoms():
					atom.transform(rotmatrix, translation)				
	current_chain+=1
	return(current_chain,ref_structure)

###################################################################################################
###################################################################################################

#Iterate of all chains in the model in order to find all residues
chain_ref='A'
chain_sample='A'

#contains a number that will be converted into a letter
current_chain=67

#Set of tuples!! RAISE ERROR IF IN A SET YOU HAVE THE SAME LETTER OR IF THIS COMBINATION iS REPEATED
#index of the chain!!

file_list=[]
name_list=[]
checked_list=set()


file_info=set()
for file in os.listdir():
	if file.endswith(".pdb"):
		x= re.split(".pdb", file)
		name_list.append(x[0][-2:])
		file_list.append(file)
		file_info.add((x[0][-2:],file))


print(file_info)
file_list.sort()
name_list.sort()


i=0
for element in name_list:
	if i == 0:
		for pair in file_info:
			if element in pair[0]:
				ref_structure=pair[1]	
	for element2 in name_list[name_list.index(element)+1:]:
		if i!=0:
			ref_structure=new_structure
		if element[0] == element2[0] and not element[0] in checked_list:
				#superimpose A with A
			for pair in file_info:
				if element2 in pair[0]:
					sample_structure=pair[1]
					chain_sample='A'
					print(element2)
					chain_ref='A'
					i+=1
					(current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)
		if element[0] == element2[1] and not element[0] in checked_list:
				#superimpose A with B
			for pair in file_info:
				if element2 in pair[0]:
					sample_structure=pair[1]
					chain_sample='B'
					print(chain_sample)
					chain_ref='A'
					i+=1
					(current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)	
		if element[1] == element2[0] and not element[1] in checked_list:
				#superimpose B with A
			for pair in file_info:
				if element2 in pair[0]:
					sample_structure=pair[1]
					chain_sample='A'
					print(chain_sample)
					chain_ref='B'
					i+=1
					(current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)
		if element[1] == element2[1] and not element[1] in checked_list:
				#superimpose B with B
			for pair in file_info:
				if element2 in pair[0]:
					sample_structure=pair[1]
					chain_sample='B'
					print(chain_sample)
					chain_ref='B'				
					i+=1
					(current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)
	checked_list.add(element[0])
	checked_list.add(element[1])

io = PDBIO()
io.set_structure(new_structure)
io.save("out.pdb")

# i=0
# file_list.sort()
# ref_structure=file_list[0]
# for file in file_list[1:]:
# 	if i!=0:
# 		ref_structure="out.pdb"
# 	i+=1
# 	sample_structure=file
# 	print(file)
# 	(current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)
# 	print(chr(current_chain))

# (current_chain,new_structure)=superimposition(ref_structure,sample_structure,chain_ref,chain_sample,current_chain)







