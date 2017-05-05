data2<- read.csv("C:/Users/benwulf/Downloads/all_uniparc_mapped_to_uniprot_and_proteomes.tdl",sep="\t")

# check for full rows
concated_data<- paste0(data2$UniParc,data2$UniProt,data2$Proteomes)
is_row_duplicated<- duplicated(concated_data)
sum(is_row_duplicated) # How many rows are duplicated.-- 8889

# check each collumn

is_UniParc_duplicated     <- duplicated(data2$UniParc)
sum(is_UniParc_duplicated)  #57719
is_UniProt_duplicated     <- duplicated(data2$UniProt)
sum(is_UniProt_duplicated)  #8889
is_Proteomes_duplicated   <- duplicated(data2$Proteomes)
sum(is_Proteomes_duplicated) #60688
is_Proteomes_duplicated_woe   <- duplicated(data2$Proteomes[data2$Proteomes!=''])
sum(is_Proteomes_duplicated_woe) #49395


# find emtpy values
sum(data2$Proteomes=='') #11294 => only collumn with empty values around 1/6 of entire data
sum(data2$UniParc=='') # 0
sum(data2$UniProt=='') #0


all.equal(is_row_duplicated,is_UniProt_duplicated) # returns TRUE => we can use UniProt as an identifier for a row


# entries in the collumn Proteomes has more than one value, they are ',' seperated

Proteomes