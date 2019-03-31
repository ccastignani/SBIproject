'''This is our main command line script to be executed'''


###################################### MODULE IMPORTS ##########################################

import os
import sys
from Bio.PDB import *
import re
from superimposer import *
from macrocomplex import *
import argparse

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

###################################### EXCEPTIONS ##########################################

class WrongInput(NameError):
	def __str__(self):
		return "ERROR: Wrong input - Only one file given"

class TooFewIterations(ValueError):
	def __init__(self,iterations):
		self.iterations=iterations
	def __str__(self):
		return "ERROR: Please try with a higher number of iterations. Minimum is 2 iterations, your input was %d." %(self.iterations)

####################################### ARGPARSER ###########################################

argument_parser = argparse.ArgumentParser(description='Build macrocomplexes')
argument_parser.add_argument('-c', '--chain_names', dest="chains",
                    action = "store",
                    default = None,
                    help = "Specify the chain names requiered to build the macrocomplex (i.e -c XABCDE)")

argument_parser.add_argument('-i', '--iterations', dest = "iterations", type=int,
                    action = "store",
                    default = None,
                    help = "Specify the number of subunits that build the macrocomplex")

argument_parser.add_argument('-coff', '--cut-off', dest = "cut_off", type=float,
                    action = "store",
                    default = 1.1,
                    help = "Specify the minimum cut-off distance for a neighbour search (used for the search of clashes) in angstroms")

argument_parser.add_argument('-o', '--out_file', dest = "out_file",
                    action = "store",
                    default = "output",
                    help = "Subunit output file name - Without extension")

options = argument_parser.parse_args()



###################################### BUILDING A STRUCTURE ##########################################

# This part creates an entire macrocomplex if all required interactions are given
# or a subunit that will be the starting point to build the eniter macrocomplex structre. 

print("\n"+"#"*100)
print("#"*100,"\n")

# Handling the error if the number of iterations asked is smaller than 2

if options.iterations != None:
	if options.iterations < 2:
		raise TooFewIterations(options.iterations)


# Initializing some variables 
name_list=[]
checked_list=set()
chain_names_list=[]
superimpose_list=set() # creating a set where already superimposed structures are added 
out_number = 1 #number of independent output structures (in case the structure contains separate complexes)



# Handling argparser arguments to create chain_names_list with all the letters of the chains 
# that we want to consider, if none are given we go through te current directory and then 
# append all the chain letters we find 

if options.chains != None: #Chains have been specified
	chain_names = options.chains
	for character in chain_names:
		if character != " ":
			chain_names_list.append(character)
else:
	for file in os.listdir():
		if file.endswith(".pdb"):
			x= re.split(".pdb", file)
			if x[0][-2:-1] not in chain_names_list:
				chain_names_list.append(x[0][-2:-1])
			if x[0][-1:] not in chain_names_list:
				chain_names_list.append(x[0][-1:])

file_info=set() #set of tuples with the chain names [0] and filename [1] ((XA,XA.pdb),...)


# creating a list of interaction chain name pairs (name_list)

for file in os.listdir():
	if file.endswith(".pdb"):
		x= re.split(".pdb", file)
		if x[0][-2:-1] in chain_names_list and x[0][-1:] in chain_names_list:
			name_list.append(x[0][-2:]) #list of chain pair names (XA,XB...)
			file_info.add((x[0][-2:],file))

name_list.sort()

print("Considering the following interaction pairs:\n"+str(name_list)+"\n")
print("="*100)

#controling the error given if there is only one input file given 
if len(file_info) == 1:
	raise WrongInput()



# going through the name_list and checking for same chain names
# then calling our superposer function with the specific input arguments

i=0 # counter for number of added structure files 
entered = None

for element in name_list:
	if i == 0: # in the first iteration. first file = reference structure file.
		current_chain1=65
		current_chain2=67
		for pair in file_info: 
			if element in pair[0]: #if element (i.e XA) in first position of tuple, retrieve file corresponding to that interacting pair
				ref_structure=pair[1]	

	# avoids to add structures that have no relevance with the complex
	if len(chain_dict)!=0 and element[0] not in checked_list and element[1] not in checked_list:			
		continue

	for element2 in name_list:
		if element != element2:
			if i!=0: 
				#if is not the first round, get the result from previous superimposer as reference structure
				ref_structure=new_structure


			#next we check each element from reference interaction pair and sample interaction pair	
			temp_str = ''.join(superimpose_list)

			if element[0] in temp_str and element[1] in temp_str and element2[0] in temp_str and element2[1] in temp_str: 
				#doesn't superimpose redundant information = when all the chains are already in the superimpose list continue
				continue

			if element[0] == element2[0] and element[0] not in checked_list: #e.g: (XA,XB) -> first chain name, matches first
				print("We are superimposeing:", element, "-", element2)
				if i!=0:
					if element in superimpose_list:				
						for pair in file_info:
							if element2 in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[0],
																										 chain_sample=element2[0],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element2)
								print("Structure",pair[1],"has been added to the model.")
								entered = True
					if element2 in superimpose_list and entered == None:	
						for pair in file_info:
							if element in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[0],
																										 chain_sample=element2[0],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element)
								print("Structure",pair[1],"has been added to the model.")
				elif i==0: 
					for pair in file_info:
						if element2 in pair[0]:	
							i+=1
							(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																									 struc_sample=pair[1],
																									 chain_ref=element[0],
																									 chain_sample=element2[0],
																									 current_chain1=current_chain1,
																									 current_chain2=current_chain2,
																									 chain_dict=chain_dict,
																									 cut_off=options.cut_off)
							superimpose_list.add(element)
							superimpose_list.add(element2)
							print("Structure",pair[1],"has been added to the model.")


			if element[0] == element2[1] and element[0] not in checked_list: #e.g: (XA,BX) -> first chain name, matches second 
				print("We are superimposing:", element, "-", element2)
				if i!=0:		
					if element in superimpose_list:				
						for pair in file_info:
							if element2 in pair[0]:	
								i+=1
								print(pair[0], pair[1])
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[0],
																										 chain_sample=element2[1],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element2)
								print("Structure",pair[1],"has been added to the model.")
								entered = True
					if element2 in superimpose_list and entered == None:	
						for pair in file_info:
							if element in pair[0]:	
								i+=1
								print(pair[0], pair[1])
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[0],
																										 chain_sample=element2[1],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element)
								print("Structure",pair[1],"has been added to the model.")						
				elif i==0:
					for pair in file_info:
						if element2 in pair[0]:	
							i+=1
							(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																									 struc_sample=pair[1],
																									 chain_ref=element[0],
																									 chain_sample=element2[1],
																									 current_chain1=current_chain1,
																									 current_chain2=current_chain2,
																									 chain_dict=chain_dict,
																									 cut_off=options.cut_off)
							superimpose_list.add(element)
							superimpose_list.add(element2)
							print("Structure",pair[1],"has been added to the model.")


			if element[1] == element2[0] and not element[1] in checked_list: #e.g: (AX,XB) -> second chain name, matches first
				print("We are superimosing:", element, "-", element2)
				if i!=0:
					if element in superimpose_list:				
						for pair in file_info:
							if element2 in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[1],
																										 chain_sample=element2[0],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element2)
								print("Structure",pair[1],"has been added to the model.")
								entered = True
					if element2 in superimpose_list and entered == None:	
						for pair in file_info:
							if element in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[1],
																										 chain_sample=element2[0],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element)
								print("Structure",pair[1],"has been added to the model.")	
				elif i==0:
					for pair in file_info:
						if element2 in pair[0]:	
							i+=1
							(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																									 struc_sample=pair[1],
																									 chain_ref=element[1],
																									 chain_sample=element2[0],
																									 current_chain1=current_chain1,
																									 current_chain2=current_chain2,
																									 chain_dict=chain_dict,
																									 cut_off=options.cut_off)
							superimpose_list.add(element)
							superimpose_list.add(element2)
							print("Structure",pair[1],"has been added to the model.")


			if element[1] == element2[1] and element[1] not in checked_list: #e.g: (AX,BX) -> second chain name, matches second
				print("We are superimposing:", element, "-", element2)
				if i!=0:
					if element in superimpose_list:				
						for pair in file_info:
							if element2 in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[1],
																										 chain_sample=element2[1],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element2)
								print("Structure",pair[1],"has been added to the model.")
								entered = True
					if element2 in superimpose_list and entered == None:	
						for pair in file_info:
							if element in pair[0]:	
								i+=1
								(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																										 struc_sample=pair[1],
																										 chain_ref=element[1],
																										 chain_sample=element2[1],
																										 current_chain1=current_chain1,
																										 current_chain2=current_chain2,
																										 chain_dict=chain_dict,
																										 cut_off=options.cut_off)
								superimpose_list.add(element)
								print("Structure",pair[1],"has been added to the model.")
				elif i==0:
					for pair in file_info:
						if element2 in pair[0]:	
							i+=1
							(current_chain1,current_chain2,chain_dict,new_structure)=superimposition(struc_ref=ref_structure,
																									 struc_sample=pair[1],
																									 chain_ref=element[1],
																									 chain_sample=element2[1],
																									 current_chain1=current_chain1,
																									 current_chain2=current_chain2,
																									 chain_dict=chain_dict,
																									 cut_off=options.cut_off)
							print("Structure",pair[1],"has been added to the model.")
							superimpose_list.add(element)
							superimpose_list.add(element2)
		entered = None
	checked_list.add(element[1])
	checked_list.add(element[0])
	
# creating the subunit output file 	
out_name= options.out_file+".cif"
io = MMCIFIO()
io.set_structure(new_structure)
io.save(out_name)

print("="*100)
print("\nStructure", options.out_file+".cif", "has been created – END OF INPUT FILES\n")
print("#"*100)


###################################### BUILDING A MACROCOMPLEX FROM SUBUNIT ##########################################

""" Here we will build a macrocomplex out of multiple subunits if asked.
	The number of iterative addition of the subunit is contolled by the interation argument -i
"""

# If the option number of iterations is given as an argument, the subunit is read from the outputed cif file and
# a list of pairs of chains with 95% or more identical sequence is built with the matcher function.
	
if options.iterations != None:
	print("\nBasic subunit built")
	print("\nCreating a macrocomplex with",str(options.iterations+1),"subunits.\n")
	print("="*100,"\n")
	(subunit,match_list) = matcher(options.out_file+".cif")

	chain_newlist=[]
	for chain in subunit[0]:
		chain_newlist.append(chain.get_id())

	# The number of the current chain is started in consistency with what has been read from the file containing the subunit
	
	current_chain1=ord(chain_newlist[-1][:1]) # First letter of last chain
	current_chain2= ord(chain_newlist[-1][1:]) # Second letter of last chain
	current_chain2+=1

	putative_structures=[] # We start a list of putative structures that will store the structures that are built from the number of iterations provided

	while len(match_list)>0:
		start=copy.deepcopy(subunit)	# The original subunit is set to the subunit from the file.
		(validation)=trial(new_structure=start,n=options.iterations,
							current_chain1=current_chain1,
							current_chain2=current_chain2,
							match_list=match_list,
							start=subunit,
							cut_off=options.cut_off)
		if validation!=False:
			new_structure=validation	# If a structure is created successfully with the required number of iterations, it is stored in the list
			putative_structures.append(new_structure)
	
	if len(putative_structures)==0: # If the number of structures created is 0, then it is recomended to lower the number of iterations
		print("============NO STRUCTURE WAS BUILT WITH SUCCESS, TRY WITH A SMALLER NUMBER OF ITERATIONS============")
	
	else:
		c=1 # Start a counter to rename correctly the file		
		# If the list of structures contains at leas one structure, they are outputet as CIF files with a given ID number
		for structure in putative_structures:
			io=MMCIFIO()
			io.set_structure(structure)
			io.save(options.out_file+"_complex_"+str(c)+".cif")
			print("Complex out_complex"+options.out_file+str(c)+".cif has been created")
			c+=1
	print("\n --- END OF PROGRAM ---")	
else:
	"\n --- END OF PROGRAM ---\n"
print("\n"+"#"*100)
print("#"*100,"\n")

