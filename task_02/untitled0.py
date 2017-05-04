# -*- coding: utf-8 -*-
"""
Created on Thu May  4 09:31:40 2017

@author: Ben
"""

## Bottom UP Version

def minsize(dictionary):
    
    mins=float('inf')
    keys=''
    for key, value in dictionary.items():
        l=len(value)
        if l<mins:
            mins=l
            keys=key
    return(keys)


def reduce(Proteom_dict,UniParc_Dict):
    required=list()
    while(len(Proteom_dict)):
        p=minsize(Proteom_dict)
        b=True
        for UniParc in Proteom_dict(p):
            if len(UniParc_Dict[UniParc])<=1:
                b=False
                break
        if(b):
            for UniParc in Proteom_dict(p):
                UniParc_Dict[UniParc].remove(p)
            
        else:
            required.append((p,Proteom_dict[p]))
        Proteom_dict.pop(p)



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
            required.append((p,unique))
            
        Proteom_dict.pop(p)
    return(required)