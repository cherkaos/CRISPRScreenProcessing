screenProcessing<-function(file,controlStart="X1_1",controlEnd="X24_2",maxsgRNA=15,minReadCount=100,zscore=FALSE){
  
#  file="APK-1-and-2-final.txt"
#  controlStart="sample35"
#  controlEnd="sample19.II"
#  maxsgRNA=15
#  minReadCount=100
#  zscore=TRUE

  data=read.csv(file,sep = "\t")
  if(zscore==TRUE){
    for(id in 3:length(data)){
      data[,id]=scale(data[,id],center=TRUE,scale=TRUE)[,1]
    }
  }

  # data=read.csv('C:/Users/sarahche/Dropbox/SARI_TILL_SCIENCE/files_for_zscore/wt_for_zscores/absolut_reads/wtall_sort_for_tgfbr.count.txt',sep = "\t")

  
  
  controlStartPos=match(controlStart,names(data))
  controlEndPos=match(controlEnd,names(data))
  ALLfoldchange=c()  
  output=data.frame()
  warning('off','all')
  for(id in controlStartPos:controlEndPos){
    foldchange=c()
    x=data[,id];
    xsorted <- data[order(x,decreasing = TRUE),c(1,2,id)] 
    xsorted0=xsorted
    #xsorted0=xsorted[xsorted[,3] != 0, ] 
    

    for(j in 1:maxsgRNA)
    {
      foldchange=c(foldchange,xsorted0[j,3]/xsorted0[j+1,3])
    }
    ALLfoldchange=c(ALLfoldchange,max(foldchange))
    posFC=match(max(foldchange),foldchange) 
    Newsorted=rbind(xsorted0[1:posFC,],Separation=c("-","-","-"),xsorted0[(posFC+1):(dim(xsorted0)[1]),])
    
    #output=c(output,Newsorted)
    if(nrow(output)==0){
      output=Newsorted
    }
    else{
      output=cbind(output,Newsorted)
    }
  }
  
  print('Hi Till, here are the FC calculated for your control clones :')
  print(ALLfoldchange)
  print('For the other clones, I will use the smalest FC as a threshold  :')
  print(min(ALLfoldchange))
  print('Which is control clone :')
  print(names(data)[match(min(ALLfoldchange),ALLfoldchange)+2])
  
  refFoldchange=min(ALLfoldchange)
  clonesToRemove=c()
  position=controlEndPos+1;
  for(i in position:length(data)){
    foldchange=c()
    x=data[,i]
    xsorted =data[order(x,decreasing = TRUE),c(1,2,i)] 
    if(!xsorted[1,3]==0){
      xsorted0=xsorted[xsorted[,3] != 0, ] 
      #To find sudden drop
      for(j in 1:maxsgRNA)
      {
        foldchange=c(foldchange,xsorted0[j,3]/xsorted0[j+1,3])
      }
      posFC=match(max(foldchange),foldchange) 
      if(max(foldchange)>refFoldchange && xsorted0[posFC,3]>minReadCount){
          
          Newsorted=rbind(xsorted0[1:posFC,],Separation=c("-","-","-"),xsorted0[(posFC+1):(dim(xsorted0)[1]),])
          
          #output=c(output,Newsorted)
          if(nrow(output)==0){
            output=Newsorted
          }
          else{
            output=cbind(output,Newsorted)
          }
          
        } 
      else {
           clonesToRemove=c(clonesToRemove,names(data)[i])
        }
    }
    else {
        clonesToRemove=c(clonesToRemove,names(data)[i])
    }
  }

  # print file adding 
  newFile=sub(".txt","_Processed.txt",file)
  write.table(output,newFile, sep="\t",row.names = FALSE,na = " ") 
}

  