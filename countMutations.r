countMutations = function(cds){
#input: a cds sequence (string)
#output: a list of three matrices with all possible Silent, Missense, and Nonsense mutation counts by mutation type and context

allmuts = c("T>G_GTT","T>G_CTT","T>G_ATT","T>G_TTT","T>G_TTG","T>G_GTG","T>G_CTG","T>G_ATG","T>G_TTC","T>G_GTC","T>G_CTC","T>G_ATC","T>G_TTA","T>G_GTA","T>G_CTA","T>G_ATA","T>C_TTT","T>C_GTT","T>C_CTT","T>C_ATT","T>C_TTG","T>C_GTG","T>C_CTG","T>C_ATG","T>C_TTC","T>C_GTC","T>C_CTC","T>C_ATC","T>C_TTA","T>C_GTA","T>C_CTA","T>C_ATA","T>A_TTT","T>A_GTT","T>A_CTT","T>A_ATT","T>A_TTG","T>A_GTG","T>A_CTG","T>A_ATG","T>A_TTC","T>A_GTC","T>A_CTC","T>A_ATC","T>A_TTA","T>A_GTA","T>A_CTA","T>A_ATA","C>T_TCT","C>T_GCT","C>T_CCT","C>T_ACT","C>T_TCG","C>T_GCG","C>T_CCG","C>T_ACG","C>T_TCC","C>T_GCC","C>T_CCC","C>T_ACC","C>T_TCA","C>T_GCA","C>T_CCA","C>T_ACA","C>G_TCT","C>G_GCT","C>G_CCT","C>G_ACT","C>G_TCG","C>G_GCG","C>G_CCG","C>G_ACG","C>G_TCC","C>G_GCC","C>G_CCC","C>G_ACC","C>G_TCA","C>G_GCA","C>G_CCA","C>G_ACA","C>A_TCT","C>A_GCT","C>A_CCT","C>A_ACT","C>A_TCG","C>A_GCG","C>A_CCG","C>A_ACG","C>A_TCC","C>A_GCC","C>A_CCC","C>A_ACC","C>A_TCA","C>A_GCA","C>A_CCA","C>A_ACA","A>C_AAC","A>C_AAG","A>C_AAT","A>C_AAA","A>C_CAA","A>C_CAC","A>C_CAG","A>C_CAT","A>C_GAA","A>C_GAC","A>C_GAG","A>C_GAT","A>C_TAA","A>C_TAC","A>C_TAG","A>C_TAT","A>G_AAA","A>G_AAC","A>G_AAG","A>G_AAT","A>G_CAA","A>G_CAC","A>G_CAG","A>G_CAT","A>G_GAA","A>G_GAC","A>G_GAG","A>G_GAT","A>G_TAA","A>G_TAC","A>G_TAG","A>G_TAT","A>T_AAA","A>T_AAC","A>T_AAG","A>T_AAT","A>T_CAA","A>T_CAC","A>T_CAG","A>T_CAT","A>T_GAA","A>T_GAC","A>T_GAG","A>T_GAT","A>T_TAA","A>T_TAC","A>T_TAG","A>T_TAT","G>A_AGA","G>A_AGC","G>A_AGG","G>A_AGT","G>A_CGA","G>A_CGC","G>A_CGG","G>A_CGT","G>A_GGA","G>A_GGC","G>A_GGG","G>A_GGT","G>A_TGA","G>A_TGC","G>A_TGG","G>A_TGT","G>C_AGA","G>C_AGC","G>C_AGG","G>C_AGT","G>C_CGA","G>C_CGC","G>C_CGG","G>C_CGT","G>C_GGA","G>C_GGC","G>C_GGG","G>C_GGT","G>C_TGA","G>C_TGC","G>C_TGG","G>C_TGT","G>T_AGA","G>T_AGC","G>T_AGG","G>T_AGT","G>T_CGA","G>T_CGC","G>T_CGG","G>T_CGT","G>T_GGA","G>T_GGC","G>T_GGG","G>T_GGT","G>T_TGA","G>T_TGC","G>T_TGG","G>T_TGT")
GENETIC_CODE = c("F","F","L","L","S","S","S","S","Y","Y","*","*","C","C","*","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D","E","E","G","G","G","G")
names(GENETIC_CODE) = c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG")

stmp  = strsplit(allmuts,split='_')
mutation_table = t(matrix(unlist(stmp),nrow=2))
colnames(mutation_table) = c('mutation','context')

contexts = unique(mutation_table[,2])

mutation_types = c('C>A','C>G','C>T','T>A','T>C','T>G','G>T','G>C','G>A','A>T','A>G','A>C')
mutation_to = substr(mutation_types,3,3)
mutation_from = substr(mutation_types,1,1)

countsSilent = matrix(0,nrow=length(mutation_types),ncol=length(contexts))
colnames(countsSilent) = contexts
rownames(countsSilent) = mutation_types
countsMissense = matrix(0,nrow=length(mutation_types),ncol=length(contexts))
colnames(countsMissense) = contexts
rownames(countsMissense) = mutation_types
countsNonsense = matrix(0,nrow=length(mutation_types),ncol=length(contexts))
colnames(countsNonsense) = contexts
rownames(countsNonsense) = mutation_types

for (iii in 1:12){
    stmp = strsplit(cds,split='')[[1]]
    if (length(stmp)%%3!=0){stop('CDS sequence length is not a multiple of 3')}
    Itmp = which(stmp==mutation_from[iii])
    Itmp1 = which(Itmp>3 & Itmp<length(stmp)-2)
    Itmp = Itmp[Itmp1]
    Itmp1 = floor((Itmp-1)/3)+1
    #Context
    contmp = paste(stmp[Itmp-1],stmp[Itmp],stmp[Itmp+1],sep='')
    codons = paste(stmp[(Itmp1-1)*3+1],stmp[(Itmp1-1)*3+2],stmp[(Itmp1-1)*3+3],sep='')
    poswithincodon = Itmp-3*(Itmp1-1)
    #Induce mutations
    codonsChanged = codons
    for (k in 1:length(codons)){
      stmp2 = codons[k]
      substr(stmp2,poswithincodon[k],poswithincodon[k]) = mutation_to[iii]
      codonsChanged[k]  =stmp2
    }
    #Now find the type of mutation
    type1 = GENETIC_CODE[match(codons,names(GENETIC_CODE))] #Before mutation
    type2 = GENETIC_CODE[match(codonsChanged,names(GENETIC_CODE))] #After mutation
    
    issilent = 1*(type1==type2)
    ismissense = 1*(type1!=type2 & type2!='*')
    isnonsense = 1*(type1!='*' & type2=='*')
    
    t = table(issilent,contmp)
    if (dim(t)[1]==1){t = rbind(t,0)}
    Itmp = match(contexts,colnames(t))
    t = t[,Itmp]
    t = as.numeric(t[2,])
    t[is.na(t)] = 0   
    countsMissense[iii,] = t
    
    t = table(ismissense,contmp)
    if (dim(t)[1]==1){t = rbind(t,0)}
    Itmp = match(contexts,colnames(t))
    t = t[,Itmp]
    t = as.numeric(t[2,])
    t[is.na(t)] = 0   
    countsMissense[iii,] = t
    
    t = table(isnonsense,contmp)
    if (dim(t)[1]==1){t = rbind(t,0)}
    Itmp = match(contexts,colnames(t))
    t = t[,Itmp]
    t = as.numeric(t[2,])
    t[is.na(t)] = 0   
    countsNonsense[iii,] = t
    
  }
  return(list(Silent=countsSilent, Missense=countsMissense, Nonsense=countsNonsense))
}


