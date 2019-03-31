from Bio.PDB import *
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import copy
from superimposer import clash_search

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def matcher(structure_cif):
	"""This function takes a file in cif format and returns a list with the pairs of chains sharing >95% of identity and the structure
	
	Input: filename of a file in cif format containing the structure that we want to use to build the macrocomplex

	Output: 
	- match_list: a list of lists containing the pairs that share enough homology to be superimposed
	- subunit: the structure parsed from the input file
	
	"""
	checkedlist=[]
	match_list=[]

	# 1. Parsing the file to get the structure

	parser = MMCIFParser()
	subunit = parser.get_structure("reference", structure_cif)

	# 2. Iterating in the chains of the subunit to do the local alignment and see how much identity they share

	for chain in subunit[0]:
		for chain2 in subunit[0]:
			if chain2.get_id() in checkedlist:  #To AVOID getting the identities in both ways e.g.: AA->AB and AB->AA
				continue
			if chain != chain2:		#For every pair of chains in the subunit, getting their aminoacid sequences
				pp=PPBuilder()
				pp2=PPBuilder()
				for pp in pp.build_peptides(chain):
					seq1 = pp.get_sequence()
				for pp2 in pp2.build_peptides(chain2):
					seq2 = pp2.get_sequence()

				# 3. Calculating the alignment and getting the score of the alignment and its length

				alignments = pairwise2.align.localxs(seq1, seq2,-1,-1)
				score= alignments[0][2]
				align_len = alignments[0][4]-alignments[0][3]

				# 4. Choosing which ones of the pairs pass the threshold and are added to the list of possible superimpositions
				
				if score == 0 or score/align_len <= 0.95: # Skipping those pairs with less than 95% identity
					continue
				else:
					match_list.append([chain.get_id(),chain2.get_id()]) # Adding to the match_list those sequences that pass the threshold
		
		checkedlist.append(chain.get_id())		
	return (subunit,match_list)

def trial(new_structure,n,current_chain1,current_chain2,match_list,start,cut_off):
	"""This function takes two structures, a number of iterations and a match_list (among others) and returns "False" or a macrocomplex

	Input:
	- new_structure: this will be the reference structure and is the original subunit
	- start: will be the sample structure in the superimposition process and it is the original subunit
	- n: number of iterations given by the user to build the macrocomplex. This is the number of subunits that form the complex -1.
	- current_chain1/2: contain the integer that will be converted to a character and used to add new chains to the structure
	- match_list: contains the list of pairs of chains that could be superimposed
	- cut_off: the minimum distance between the atoms of a new chain and the existing structure that will allow the chain to be added


	A new trial is started while the match_list still contains pairs of chains that could be superimposed. A trial ends when an attempt
	to superimpose two chains fail, those two chains need to be removed from the match_list, to not be tried again. A trial also ends when 
	the number of superimpositions reaches the number of iterations indicated by the user. Also in this case the matc_list is actualized.


	Output:
	- FALSE will be returned when more superimpositions cannot be performed with a given pair of chains and the number of iterations was not reached
	- NEW_STRUCTURE will be returned if a pair of chains reaches the number of superimpositions required

	"""
	
	# 1. Deepcopy of the original subunit to restart the coordinates of the sample subunit 
	#    Start a chain dictionary that contains the original name of a chain in the subunit as key, 
	#    and all the new names that this chain has acquired in order to be added to the structure as values. 
	#    The superimposed chain that doesn't need to be added to the structure (is reseted in each trial)
	#	 Start the counter of performed iterations as j=1. 
	
	original=copy.deepcopy(start)
	chain_dict={} 													
	skip_chain=None					
	j=1  							
	
	# 2. Calling the function jump that will produce two possible outputs, False when a superimposition failed, or a new structure when it worked
	# 	 Taking different values from the function that are to keep track of different values that need to be followed
	#    Taking Validation that contains either False or a structure
	
	while j<=n:
		(validation,current_chain1,current_chain2,skip_chain,chain_dict)=jump(new_structure=new_structure,
																			j=j,current_chain1=current_chain1,
																			current_chain2=current_chain2,
																			skip_chain=skip_chain,
																			chain_dict=chain_dict,n=n,
																			original=original,match_list=match_list,
																			cut_off=cut_off)
		if validation==False:	# If validation is False, then the superimposition didn't work and a new trial has to start
			return(False)
		
		else:					# If validation is not false, then the new structure is saved and j increases 1						
			j+=1
			new_structure=validation

	print("\n"+"="*35,"STRUCTURE CREATED SUCCESSFULLY","="*35,"\n")
	return(new_structure)  # If the number of iterations is reached, the structure is returned

def jump(new_structure,original,j,current_chain1,current_chain2,skip_chain,chain_dict,n,match_list,cut_off):
	"""This function takes two structures, a match_list (among others) and returns False or a new structure (among others)

	Input:
	- new_structure: the reference structure (in the first iteration it is the original one, in further iterations it is the amplified one)
	- original: will be the sample structure in the superimposition process and it is the original subunit
	- n: number of iterations given by the user to build the macrocomplex. This is the number of subunits that form the complex -1.
	- current_chain1/2: contain the integer that will be converted to a character and used to add new chains to the structure
	- match_list: contains the list of pairs of chains that could be superimposed
	- cut_off: the minimum distance between the atoms of a new chain and the existing structure that will allow the chain to be added
	- skip_chain: in the first iteration it is set to the reference chain in the superimposition, then it is maintained
	- chain_dict: dictionary containing the original name of a chain in the structure as key and all the new names it acquires as values
	- j: number of the current iteration

	Every JUMP is an attempt to superpose two chains that are in the first item of the match_list. When those chain matching the list are found
	the superimposition is tried and from its return value this function will perform different operations. 
	If the return is False, meaning that the superimposition didn't work, the match is removed from the list and False is returned. 
	If the return value is not false, then the new structure is returned and the match list is actualized to the new names of the matching chains

	Output:
	- FALSE will be returned when more superimpositions cannot be performed with a given pair of chains and the number of iterations was not reached
	- NEW_STRUCTURE will be returned if a pair of chains reaches the number of superimpositions required

	Both FALSE and the NEW_STRUCTURE are returned with other values that the functions need to keep track of
	
	"""

	# 1. Iterating through the chains in the original (sample) and the new_structure (reference) to find those chains taht match the first item on
	# 	 the match list. When they are found, the SUPERIMPOSE function is called.  

	for chain in original[0]:
		for chain2 in new_structure[0]:
			if chain != chain2:
				match=match_list[0]
				if chain.get_id()==match[0] and chain2.get_id() == match[1]:  # Checking if the chains match the first pair on the list
					if j==1:
						skip_chain=chain2.get_id()		# In the first iteration we set the skip chain to the one in the reference structure

					(validation,current_chain1,current_chain2,skip_chain,chain_dict)=superimpose(ref_structure=new_structure,ref_chain=chain2,n=n,
																									sample_chain=chain,skip=skip_chain,
																									current_chain1=current_chain1,
																									current_chain2=current_chain2,
																									chain_dict=chain_dict,j=j,
																									original=original,cut_off=cut_off)
					if validation==False:	# If the superimposition didn't work, we remove the pair from the list and return False
						match_list.pop(0)
						return (False,current_chain1,current_chain2,skip_chain,chain_dict)
					
					else:

						# In the first iteration, the match_list is actualized with the name of the reference chain in the left position
						# (this will be found on the sample structure in further iterations) and in the right position the new name of the
						# chain that matched the reference one (which is only found in the structure that is growing)
						
						if j==1:			
							match_list[0]=[chain2.get_id(),chain_dict[chain.get_id()][j-1]]
						
						# If the number of iterations is reached, the pair is removed from the list
						
						elif j==n:
							match_list.pop(0)
						
						# If it is not the first nor the last iteration, the value on the left is the chain in the sample structure 
						# (that will be the same on further iterations, and found in the sample structure) and the value on the right
						# the new name that the reference chain has taken
						
						else:
							for key,value in chain_dict.items():
								if chain2.get_id() in value:
									match_list[0]=[chain.get_id(),chain_dict[key][j-1]]

						# After actualizing the match_list, the new_structure is returned
						new_structure=validation  
						return (new_structure,current_chain1,current_chain2,skip_chain,chain_dict)

def superimpose(ref_structure,original,n,ref_chain,sample_chain,skip,current_chain1,current_chain2,chain_dict,j,cut_off):
	"""This function takes reference and sample structures and chains (among others) and returns false or a new structure (among others)

	Input:
	- ref_structure: the reference structure (original subunit or amplified one depending on the iteration number)
	- original: will be the sample structure in the superimposition process and it is the original subunit
	- n: number of iterations given by the user to build the macrocomplex. This is the number of subunits that form the complex -1.
	- ref_chain: chain from the reference structure to be superimposed
	- sample_chain: chain from the sample structure to be superimposed
	- skip: chain that needs to be skipped when adding chains to the structure with transformed atoms (in the first iteration it is a chain belonging
	to the original subunit)
	- current_chain1/2: contain the integer that will be converted to a character and used to add new chains to the structure
	- match_list: contains the list of pairs of chains that could be superimposed
	- chain_dict: dictionary containing the original name of a chain in the structure as key and all the new names it acquires as values
	- j: number of the current iteration
	- cut_off: the minimum distance between the atoms of a new chain and the existing structure that will allow the chain to be added

	The superimpose tries to superpose the chains that are given by the jump function, it will try to add the chains from the sample structure
	to the one that is being amplified. Then it will check for clashes calling a clash_check function. If the return value is True, it detaches
	the new chain that had been added to the structure and returns False and other variables to the JUMP function. If the return value is False
	for all the chains that need to be added, then the function will return a deepcopy of the structure amplified among other variables.
	
	Output:
	- FALSE will be returned when more superimpositions cannot be performed with a given pair of chains and the number of iterations was not reached
	- NEW_STRUCTURE will be returned if a pair of chains reaches the number of superimpositions required

	Both FALSE and the NEW_STRUCTURE are returned with other values that the functions need to keep track of

	"""

	print("Superimposing",sample_chain.get_id(),"===",ref_chain.get_id())
	
	# Initializing two lists, one for the set of reference atoms and one for the set of sample atoms, sample atoms are the ones to which
	# we will apply the superimposition (the ones that move)
	ref_atoms=[]
	sample_atoms=[]
	
	# 1. First it is checked if the length of the polypeptide sequence of the two chains is equal or not. If it is not equal, we will take 
	#    only the atoms from the shared residues. If lengths are the same, we take all atoms.

	pp=PPBuilder()
	pp2=PPBuilder()
	for pp in pp.build_peptides(ref_chain):
		seq_ref=pp.get_sequence()
	for pp2 in pp2.build_peptides(sample_chain):
		seq_sample=pp2.get_sequence()
	

	# Building the proper conditions in order to move in the sequence that is longer until we find matching residues to take the atoms
	# We move to the next residue in the longer sequence if the residue name does not match. When both chains share a residue, its CA atoms are taken

	if len(seq_ref)>len(seq_sample):
		for samp_res in sample_chain.get_residues():
			for ref_res in ref_chain.get_residues():
				if ref_res.get_id()[1] != samp_res.get_id()[1]:
					continue
				else:
					ref_atoms.append(ref_res['CA'])
					sample_atoms.append(samp_res['CA'])

	elif len(seq_sample)>len(seq_ref):
		for ref_res in ref_chain.get_residues():
			for samp_res in sample_chain.get_residues():
				if ref_res.get_id()[1] != samp_res.get_id()[1]:
					continue
				else:
					ref_atoms.append(ref_res['CA'])
					sample_atoms.append(samp_res['CA'])
	
	# If both sequences share length, we take all the residues that don't have the "H_" tags in the first position of the ID, to remove the HETATOMs

	else:
		for residue in ref_chain:
			if residue.get_id()[0] == " ":
				ref_atoms.append(residue['CA'])

		for residue in sample_chain:
			if residue.get_id()[0] == " ":
				sample_atoms.append(residue['CA'])	

	# 2. Defining a Superimposer object, we set the fixed (ref_atoms) and moving (sample_atoms), and apply the superimposition to the sample atoms.
	#    Then the rotation and translation matrices are stored in order to perform the transformation of the atoms that are added to the structure
	
	super_imposer = Superimposer()
	super_imposer.set_atoms(ref_atoms, sample_atoms)
	super_imposer.apply(sample_atoms)

	rotmatrix = super_imposer.rotran[0]
	translation = super_imposer.rotran[1]
	
	# 3. In the first iteration the skip chain was set to be the reference one (as it is the one from the sample that will not be added to the structure,
	# 	 it doesn't need to be added because its information is already in the subunit to which things are being added). e.g.: if A and B are being 
	# 	 superposed, note that the subunits are both ABC and ABC, then when superposing A on top of B, the information is repeated, so we don't need to add B. 
	# 	 In the last iterations the addition of whole subunit in those structures that are cylindric would cause clashes, we need to also avoid
	# 	 adding the sample chain (in the example that would be A) to avoid clashes and be able to buil the complete structure.
	
	if j==n:
		for key,value in chain_dict.items():
			if ref_chain.get_id() in value:
				skip2=key
	# If it is not the las iteration, we skip just one chain
	else:
		skip2=None
	
	# 4. We iterate on the original structure and for those chains that are not the ones to skip, we add a new chain taking the current chain names
	#    to be followed. Then we take the residues from that chain in the original subunit and (if they are not HETATOMs) we add them to the chain that
	#    we have added to the reference structure (the one to be amplified). Then, the atoms are transformed according to the rotation and translation
	#    matrices that were obtained. 

	for chain in original[0]:
		if chain.get_id() != skip and chain.get_id() != skip2:
			new_chain=Chain.Chain(chr(current_chain1)+chr(current_chain2))	
			ref_structure[0].add(new_chain)		# Add the chain to the reference structure
			for residue in chain.get_residues():
				if residue.get_id()[0] == " ":
					res=residue.get_id()[1]
					ref_structure[0][chr(current_chain1)+chr(current_chain2)].add(residue)		# Add the residues to the reference structure
					for atom in ref_structure[0][chr(current_chain1)+chr(current_chain2)][res].get_atoms(): # Get the atoms added to reference and transform
						atom.transform(rotmatrix, translation)

			# 5. A check for clashes is performed. If the return value of the clash_check function is True (Clashes were found between)
			#    the new chain and the chains that were already in the structure), then that chain is detached from the structure and 
			#    FALSE is returned among other variables to the JUMP function.
			#    If the return value of the clash_check is False (meaning that clashes were not found), the new chain name is added to the 
			#    dictionary (appended to the list of values that an original chain has acquired). And another chain tries to be added
			#    until there are no chains in the original strucure (different from the ones skipped on purpose) that have not been added.

			clash_check=clash_search(ref_structure,chr(current_chain1)+chr(current_chain2),cut_off=cut_off)
			if clash_check == True:
				ref_structure[0].detach_child(new_chain.get_id()) # New chain has clashes with existing ones, detaching the new chain
				return (False,current_chain1,current_chain2,skip,chain_dict) #Chain is not added, False is returned
			elif clash_check == False:
				chain_dict.setdefault(chain.get_id(),[]) # In case it is the first iteration, to avoid KeyErrors
				chain_dict[chain.get_id()].append(chr(current_chain1)+chr(current_chain2))
				if current_chain2 == 90:
					current_chain2 = 65
					current_chain1+=1
				else:
					current_chain2+=1				
	# 6. At this point, all the chains that needed to be added had been added, a deepcopy of the reference structure amplified is returned to the JUMP
	#    function, as the superimposition worked propperly

	new_structure=copy.deepcopy(ref_structure)
	return (new_structure,current_chain1,current_chain2,skip,chain_dict)

		
