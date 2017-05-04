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
        
        # create uniProt:uniPrac dict
        if line[0] in uniProt_dict:
            uniProt_dict[line[0]].update(set([line[1]]))
        else:
            uniProt_dict[line[0]] = set([line[1]])

        # protemoe exist
        if len(line) == 3:
            #extract proteoms
            proteome={x.split(':')[0] for x in line[2].split(', ')}
            # create uniPrac:proteome dict
            if line[0] in uniPrac_dict:
                uniPrac_dict[line[0]].update(proteome)     
            else:
                uniPrac_dict[line[0]] = proteome
           
            # create proteome:uniPrac dict

            for prot in proteome:
                if prot in proteome_dict:
                    proteome_dict[prot].update(set([line[0]]))
                else:
                    proteome_dict[prot] = set([line[0]])
 
    f.close()
#
## sort proteomes by number of uniPrac in proteome
#for k in sorted(proteome_dict, key=lambda k: len(proteome_dict[k])):
#        proteome_sorted_list.append(proteome_dict[k])
#
#

## Top Down Version

'''Better Version of it '''
def maxsize(dictionary,covered):
    
    maxi=-1
    keys=''
    for key, value in dictionary.items():
        l=len(value.difference(covered))
        if (l>maxi):
            maxi=l
            keys=key
    return(keys)            
            
                
def reduce_max(Proteom_dict,UniParc_Dict):
    required=list()
    covered_UniPark=set()
    while(len(Proteom_dict)):
        p=maxsize(Proteom_dict,covered_UniPark)
        unique=Proteom_dict[p].difference(covered_UniPark)
        if(len(unique)):
            covered_UniPark.update(unique)
            required.append((p,unique))
            
        Proteom_dict.pop(p)
    print("%i proteoms given, reduced to %i, %.2f%% of input"%(len(proteome_dict),len(required),(len(required)*100/len(proteome_dict))))
    return(required)

zz=reduce_max(proteome_dict.copy(),uniPrac_dict.copy())

def proteome2uniprot(proteomeList,uniProt_dict):
    ret=list()
    for elem in proteomeList:
        uprotset=set()
        for upark in elem[1]:
            uprotset.update(uniProt_dict[upark])
        ret.append((elem[0],list(uprotset)))
    #ret.sort(key=lambda tup: len(tup[1]))
    return(ret)
            
