from Bio.PDB import *
import re
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# initializing variables
# dictionary of chain names of the subunit as keys and the new names of the added chains in our structure
chain_dict={} 

def clash_search(ref_structure, new_chain, cut_off):
	''' 
	This function searches for the neighbor atoms within a certain distance of a reference atom.
	It will go through the structure chain by chain and check if the new chain we want to add 
	finds neighboring atoms within the given cut-off.
	
	Input:
	- ref_structure: the reference structure that the sample structure will be added to
	- new_chain: Name of the new chain to be used as reference for the neighbor search
	- cut_off: the minimum distance between the atoms of the new_chain and the existing structure 
	that will allow the chain to be added

	Output:
	- if it finds neighbours within the cut-off it returns: TRUE
	- if not it returns: FALSE

	'''

	print("We are looking for clashes...")
	for chain in ref_structure[0]:
		clashes=[]
		if chain.get_id() != new_chain:
			atom_list=Selection.unfold_entities(ref_structure[0][chain.get_id()], 'A')
			ns = NeighborSearch(atom_list)
			for atom in ref_structure[0][new_chain].get_atoms():
				center = atom.get_coord()
				neighbors = ns.search(center, cut_off, level='C') # 1.1 for distance in angstrom
				for element in neighbors:
					if element != ref_structure[0][new_chain]:
						clashes.append(element)
			if len(clashes) != 0:
				return True
			else:
				return False	

def superimposition(struc_ref,struc_sample,chain_ref,chain_sample,current_chain1,current_chain2,chain_dict, cut_off):
	'''
	This function calculates and applies a superposition matrix of the sample structure onto the reference structure.
	Further it returns a new structure with the sample structure added to the reference structure. 

	Input:
	- struc_ref: the reference structure that the sample structure will be added to
	- struc_sample: will be the sample structure in the superimposition process and will be added to the reference structure
	- chain_ref: chain from the reference structure to be superimposed
	- chain_sample: chain from the sample structure to be superimposed
	- current_chain1/2: contain the integer that will be converted to a character and used to add new chains to the structure
	- chain_dict: dictionary containing the original name of a chain in the structure as key and all the new names it acquires as values
	- cut_off: the minimum distance between the atoms of a new chain and the existing structure that will allow the chain to be added

	
	Output:
	- NEW_STRUCTURE (with sample chains added) and the updated chain_dict dictionary 
	'''

	parser = PDBParser()

	#the first time, we give a pdb as reference file
	if type(struc_ref) is str: 
		# getting the structures
		ref_structure = parser.get_structure("reference", struc_ref)
		new_structure = Structure.Structure('test')
		new_model= Model.Model(0)
		new_structure.add(new_model)
		new_chain = Chain.Chain("AA")
		new_structure[0].add(new_chain)
		new_chain = Chain.Chain("AB")
		new_structure[0].add(new_chain)

		# defining variables
		ref_chain_list=[] #list of chains of reference file
		original_chains=set()
		chain_dict={} # we have to start a new structure in order to keep track of what we put in the new structure
		
		# renaming the chains of the reference file 
		x= re.split(".pdb", struc_ref)
		x=x[0][-2:]
		chain_dict[x[0]]='AA'	
		chain_dict[x[1]]='AB'

		for chain in ref_structure[0]:
			ref_chain_list.append(chain.get_id()) 
		for residue in ref_structure[0][ref_chain_list[0]].get_residues():
			new_structure[0]['AA'].add(residue)
		for residue in ref_structure[0][ref_chain_list[1]].get_residues():
			new_structure[0]['AB'].add(residue)			

	# after the first iteration we pass a structure instance
	else: 
		new_structure = struc_ref

	sample_structure = parser.get_structure("sample", struc_sample)

	ref_model    = new_structure[0]
	sample_model = sample_structure[0]

	sample_chain_list=[]
	for chain in sample_structure[0]:
		sample_chain_list.append(chain.get_id()) #list of chains of reference file

	# filechain_strucchain is a dictionary that saves the name relations from the reference pdb to the structure.
	# The keys are the original names in the pdb and the correspondig value the name in the newly created structure.
	# e.g. XB.pdb will correspond to  X->AA and B->AB in the dictionary. 
	filechain_strucchain={}
	y= re.split(".pdb", struc_sample)
	y=y[0][-2:]
	filechain_strucchain[y[0]]=sample_chain_list[0]
	filechain_strucchain[y[1]]=sample_chain_list[1]
	
	# getting the CA atoms of the reference and the sample structure for 
	# calculation of the rotran matrix. 
	ref_atoms = []
	sample_atoms = []
	chain_number=0
	for ref_chain in ref_model:
		if ref_chain.get_id() == chain_dict[chain_sample]:
			#Iterate of all residues in each model in order to find proper atoms
			for ref_res in ref_chain:
				if ref_res.get_id()[0] == " ": #doesn't read heteroatoms
					ref_atoms.append(ref_res['CA'])	
	for sample_chain in sample_model:
		if sample_chain.get_id() == filechain_strucchain[chain_sample]:
			for sample_res in sample_chain:
				if sample_res.get_id()[0] == " ": #doesn't read heteroatoms
					sample_atoms.append(sample_res['CA'])					

	super_imposer = Superimposer()
	super_imposer.set_atoms(ref_atoms, sample_atoms)
	super_imposer.apply(sample_model.get_atoms())
	rotmatrix = super_imposer.rotran[0]
	translation = super_imposer.rotran[1]

	# Create new chain for the sueprimposition result
	new_chain=Chain.Chain(chr(current_chain1)+chr(current_chain2))
	new_structure[0].add(new_chain)


	# Applying the rotran matrix to every atom of the sample model
	# and checking clashes. If we find clashes we remove this
	# chain that causes the clashes. 
	for chain in sample_model:	 
		if chain.get_id()!=filechain_strucchain[chain_sample]:	
			for residue in chain.get_residues():
				if residue.get_id()[0] == " ":
					res=residue.get_id()[1]
					new_structure[0][chr(current_chain1)+chr(current_chain2)].add(residue)
					for atom in new_structure[0][chr(current_chain1)+chr(current_chain2)][res].get_atoms():
						atom.transform(rotmatrix, translation)	
			clash_check=clash_search(new_structure,chr(current_chain1)+chr(current_chain2), cut_off)
			if clash_check == True:
				new_structure[0].detach_child(new_chain.get_id())
				print("!!! WE FOUND CLASHES !!!")
				print(new_chain.get_id(), "has been detached")

			if clash_check == False:
				print("NO clashes found")
				if y[0] != chain_sample and y[0] not in chain_dict.keys():
					chain_dict[y[0]]=chr(current_chain1)+chr(current_chain2)
					if current_chain2 == 90:
						current_chain2 = 65
						current_chain1+=1
					else:
						current_chain2+=1	
				elif y[1] != chain_sample and y[1] not in chain_dict.keys():
					chain_dict[y[1]]=chr(current_chain1)+chr(current_chain2)
					if current_chain2 == 90:
						current_chain2 = 65
						current_chain1+=1
					else:
						current_chain2+=1
				else:
					new_structure[0].detach_child(new_chain.get_id())

	return(current_chain1,current_chain2,chain_dict,new_structure)





