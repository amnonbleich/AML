# Machine Learning Assignment 2 - Lee Hong, Amnon Bleich, Ben Wulf

import json
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
## get the Proteom with the most included UniParc ID's
def maxsize(dictionary,covered):
    
    maxi = -1
    keys = ''
    for key, value in dictionary.items():
        l = len(value.difference(covered))
        if (l>maxi):
            maxi = l
            keys = key
    return(keys)            
 

## creates a list of Tuples (Proteom, [uniParc IDs])                        
def reduce_max(local_proteome_dict):
    required = list()
    covered_uniPrac = set()
    while(len(local_proteome_dict)):
        p = maxsize(local_proteome_dict,covered_uniPrac)
        unique = local_proteome_dict[p].difference(covered_uniPrac)
        if(len(unique)):    # if there is new information
            covered_uniPrac.update(unique) # add this info to covered uniparc
            required.append((p,unique))     # and add the proteom  to the list of required proteomes
            
        local_proteome_dict.pop(p)      # remove the proteome from the initial set
    print("%i proteomes given, reduced to %i, %.2f%% of input"%(len(proteome_dict),len(required),(len(required)*100/len(proteome_dict))))
    return(required)



zz=reduce_max(proteome_dict.copy()) 

## Remap UniParc identifier to UniProt and Replace them
def proteome2uniprot(proteomeList,uniPrac_dict):
    ret = list()
    for elem in proteomeList:       # for each proteom
        uProt_set = set()               # create a new set
        for uPrac in elem[1]:               # translate uniparc to uniprot
            uProt_set.update(uniPrac_dict[uPrac])
        ret.append((elem[0],list(uProt_set)))## Create a new list with UniProt IDs instead of unipark
    # sort proteomes by number of uniProt associated
    ret = sorted(ret , key=lambda x: len(x[1]), reverse=True) # Resort, could be that the number of included UniProt changes the relation
    return(ret)

proteome_sorted_list = proteome2uniprot(zz, uniPrac_dict)


# Write solution file
with open('./data.js', 'w') as outfile:
    json.dump(proteome_sorted_list, outfile)

#a=[len(x[1]) for x in proteome_sorted_list]    
#plt.figure(num=None, figsize=(15, 10), dpi=80, facecolor='w', edgecolor='k')
#plt.plot(a)
#plt.ylabel('Number of UniProt ID\'s')
#plt.savefig('foo.png')



#Sainity check
alle=list()
b=set()
for a in proteome_sorted_list:
    alle +=a[1]
    b.update(set(a[1]))
    
    
print(len(alle))#54326
print(len(b))#54326 -> they are unique

## 54326 unique uniprot in our mapping but 54712 in original 386 not mapped. What happend?
## It seems that there are 386 UniProt which have no reference to a Proteom
## Maybe we can create a list with them. but we think it wasn't asked