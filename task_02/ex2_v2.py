FILE_PATH = './all_uniparc_mapped_to_uniprot_and_proteomes.tdl'

# dictionary UniPrac: tuple
uniPrac2Proteome_dict = {}
uniPrac2uniProt_dict = {}
uniProt2uniParc_dict = {}
proteome2uniParc_dict = {}

uniProt2Proteome_dict={}
Proteom2uniProt_dict={}



with open(FILE_PATH) as f:
    next(f) # skipe first line
    for line in f:
        # data processing
        line = line.rstrip() # remove special characters 
        line = line.split('\t')
        
        # create uniParc:uniProt dict
        if line[0] in uniPrac2uniProt_dict:
            uniPrac2uniProt_dict[line[0]].update(set([line[1]]))
        else:
            uniPrac2uniProt_dict[line[0]] = set([line[1]])
            
            # create uniParc:uniProt dict
        if line[1] in uniProt2uniParc_dict:
            uniProt2uniParc_dict[line[1]].update(set([line[0]]))
        else:
            uniProt2uniParc_dict[line[1]] = set([line[0]])

        # protemoe exist
        if len(line) == 3:
            #extract proteoms
            proteome={x.split(':')[0] for x in line[2].split(', ')}
            # create uniPrac:proteome dict
            if line[0] in uniPrac2Proteome_dict:
                uniPrac2Proteome_dict[line[0]].update(proteome)     
            else:
                uniPrac2Proteome_dict[line[0]] = proteome
           
            # create proteome:uniPrac dict

            for prot in proteome:
                if prot in proteome2uniParc_dict:
                    proteome2uniParc_dict[prot].update(set([line[0]]))
                else:
                    proteome2uniParc_dict[prot] = set([line[0]])
 
    f.close()
    
    
'''Create uniProt2Proteome and Proteom2uniProt'''

for proteom in proteome2uniParc_dict:
    uprotset=set()
    for uniParc in proteome2uniParc_dict[proteom]:
        uprotset.update(uniPrac2uniProt_dict[uniParc])
    Proteom2uniProt_dict[proteom]=uprotset

for uprot in uniProt2uniParc_dict:
    proteomset=set()
    for upark in uniProt2uniParc_dict[uprot]:
        proteomset.update(uniPrac2Proteome_dict[upark])
          
    uniProt2Proteome_dict[uprot]=proteomset


    
    
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
    print("%i proteoms given, reduced to %i, %.2f%% of input"%(len(Proteome_dict),len(required),(len(required)*100/len(Proteome_dict))))
    return(required)

zz=reduce_max(Proteom2uniProt_dict.copy(),uniProt2Proteome_dict.copy())


