#Simulate Random Walk Algorithm
#Input parameters for restart rate and threshold value
args=(commandArgs(TRUE))
r=0.7
threshold=0.00001

#setwd('write working directory')
strTargetDicFile= './processed_data/rwr_gene_dictionary.txt'
strNetworkFile='./processed_data/rwr_relation_dictionary.txt'
strTargetInteractionFile= './processed_data/rwr_drug_target_dictionary.txt'

DBName=strsplit(strNetworkFile, ".txt")[[1]]
TargetDic <- read.delim(strTargetDicFile, header=F)
Network_size=dim(TargetDic)[1]

w=mat.or.vec(Network_size, Network_size)

Model <- read.delim(strNetworkFile, header=F)
x=Model[1]
y=Model[2]
RelType=Model[3]
EdgeWeight=Model[4]

N=dim(Model[2])[1]


for(i in 1:N){    
  x[i,]
  y[i,]
  if (RelType[i,]=="->"){
    w[x[i,],y[i,]]=as.numeric(EdgeWeight[i,])
  }
}
  

someenv<-new.env()
ModelTarget <- read.delim(strTargetInteractionFile, header=F)
N_Target=dim(ModelTarget[2])[1]
x_Target=ModelTarget[1]
y_Target=ModelTarget[2]
for(i in 1:N_Target) {    
  if (!is.na(y_Target[i,])){    
    someenv[[toString(x_Target[i,])]]=append(someenv[[toString(x_Target[i,])]],y_Target[i,])
  }  
}

keys=ls(someenv)
for(i in 1:length(keys)) {
  # print (toString(100.0*i/length(keys)))
  print (keys[i])
  if (grepl('/',keys[i])) {keys[i]=gsub('/','',keys[i])}
  foName = paste0('./random_walk_result/',keys[i],'.txt')
  p0=t(mat.or.vec(1,Network_size))  
  for (j in 1:length(someenv[[keys[i]]])){    
    p0[someenv[[keys[i]]][j]]=30
  }  
  
  p=p0     
  MatrixVariationValue=1  
  while (MatrixVariationValue>threshold){  
    before_p=p
    p=((1-r)*t(w))%*%p+r*p0
    after_p=p
    MatrixVariationValue=sum(abs(before_p-after_p))
  }  
  
  rownames(after_p) <- TargetDic$V1
  print (length(which(after_p!=0.0)))
  write.table(after_p,file=foName, row.names = TRUE, col.names=FALSE,sep="\t", quote=FALSE,eol = "\n")
}
