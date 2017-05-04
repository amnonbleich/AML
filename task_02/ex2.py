FILE_PATH = './all_uniparc_mapped_to_uniprot_and_proteomes.tdl'

# dictionary UniPrac:uniProt
uniPrac_dict = {}
# dictionary proteome:uniPrac
proteome_dict = {}

with open(FILE_PATH) as f:
    next(f) # skip first line
    for line in f:
        # data processing
        line = line.rstrip() # remove special characters 
        line = line.split('\t')
        
        # create uniPrac:uniProt dict
        if line[0] in uniPrac_dict:
            uniPrac_dict[line[0]].update(set([line[1]]))
        else:
            uniPrac_dict[line[0]] = set([line[1]])

        # protemoe exist
        if len(line) == 3:
            #extract proteoms
            proteome={x.split(':')[0] for x in line[2].split(', ')}
            # create proteome:uniPrac dict
            for prot in proteome:
                if prot in proteome_dict:
                    proteome_dict[prot].update(set([line[0]]))
                else:
                    proteome_dict[prot] = set([line[0]])
 
    f.close()


## Top Down Version

def maxsize(dictionary,covered):
    
    maxi = -1
    keys = ''
    for key, value in dictionary.items():
        l = len(value.difference(covered))
        if (l>maxi):
            maxi = l
            keys = key
    return(keys)            
                         
def reduce_max(local_proteome_dict):
    required = list()
    covered_uniPrac = set()
    while(len(local_proteome_dict)):
        p = maxsize(local_proteome_dict,covered_uniPrac)
        unique = local_proteome_dict[p].difference(covered_uniPrac)
        if(len(unique)):
            covered_uniPrac.update(unique)
            required.append((p,unique))
            
        local_proteome_dict.pop(p)
    print("%i proteomes given, reduced to %i, %.2f%% of input"%(len(proteome_dict),len(required),(len(required)*100/len(proteome_dict))))
    return(required)

zz=reduce_max(proteome_dict.copy())

def proteome2uniprot(proteomeList,uniPrac_dict):
    ret = list()
    for elem in proteomeList:
        uProt_set = set()
        for uPrac in elem[1]:
            uProt_set.update(uniPrac_dict[uPrac])
        ret.append((elem[0],list(uProt_set)))
    # sort proteomes by number of uniProt associated
    ret = sorted(ret , key=lambda x: len(x[1]), reverse=True)
    return(ret)

proteome_sorted_list = proteome2uniprot(zz, uniPrac_dict)