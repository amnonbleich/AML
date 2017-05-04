FILE_PATH = './all_uniparc_mapped_to_uniprot_and_proteomes.tdl'

# dictionary UniPrac: tuple
uniPrac_dict = {}
uniProt_dict = {}
proteome_dict = {}
proteome_sorted_list = []

with open(FILE_PATH) as f:
		next(f) # skipe first line
		for line in f:
			# data processing
			line = line.rstrip() # remove special characters 
			line = line.split('\t')
			if len(line) == 3:
				line[2] = line[2].split(':')[0] # remove proteome description
			
			# create uniProt:uniPrac dict
			if line[1] in uniProt_dict:
				uniProt_dict[line[1]].update(set([line[0]]))
			else:
				uniProt_dict[line[1]] = set([line[0]])

			# protemoe exist
			if len(line) == 3:
				# create uniPrac:proteome dict
				if line[0] in uniPrac_dict:
					uniPrac_dict[line[0]].update(set([line[2]]))
				else:
					uniPrac_dict[line[0]] = set([line[2]])

				# create proteome:uniPrac dict
				if line[2] in proteome_dict:
					proteome_dict[line[2]].update(set([line[0]]))
				else:
					proteome_dict[line[2]] = set([line[0]])
		f.close()

# sort proteomes by number of uniPrac in proteome
for k in sorted(proteome_dict, key=lambda k: len(proteome_dict[k])):
        proteome_sorted_list.append(proteome_dict[k])


