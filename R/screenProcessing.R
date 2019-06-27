#' screenProcessing function
#'
#' @param inputfile csv or txt file (comma separated) where each column represent a clones and and each row.
#' @param controlStart ..
#' @param controlEnd ..
#' @param resultfile ..
#' @param maxsgRNA ..
#' @param minReadCount ..
#' @param zscore ..
#' @param orderOutput 
#' @param shortOutput 
#'
#' @return resultfile
#' @export screenProcessing
#'
#' @examples
screenProcessing<-function(inputfile,controlStart,controlEnd,resultfile,maxsgRNA,minReadCount,zscore,orderOutput=TRUE,shortOutput=TRUE){
  
  options(warn=-1)
  if(!file.exists(inputfile)) stop('No such file "', inputfile,'"')
  
  data=utils::read.csv(inputfile,sep = "\t", check.name=FALSE)
  dataNotZscored=data
  if(zscore==TRUE){
    for(id in 3:length(data)){
      data[,id]=scale(data[,id],center=TRUE,scale=TRUE)[,1]
    }
  }
  
  controlStartPos=c();
  controlEndPos=c();
  controlStartPos=match(controlStart,names(data))
  controlEndPos=match(controlEnd,names(data))
  
  if(is.na(controlStartPos)) stop('The sample"', controlStart,'" is not part of samples "', names(data),'"')
  if(is.na(controlEndPos)) stop('The sample"', controlEndPos,'" is not part of samples "', names(data),'"')
  
  ALLfoldchange=c()
  output=data.frame()
  keepOrder=c()
  for(id in controlStartPos:controlEndPos){
    foldchange=c()
    x=data[,id];
    xsorted <- data[order(x,decreasing = TRUE),c(1,2,id)]
    xsorted0=xsorted;
    #xsorted0=xsorted[xsorted[,3] != 0, ]
    
    
    for(j in 1:maxsgRNA)
    {
      foldchange=c(foldchange,xsorted0[j,3]/xsorted0[j+1,3])
    }
    ALLfoldchange=c(ALLfoldchange,max(foldchange))
    posFC=match(max(foldchange),foldchange)
    
    if(!shortOutput){
      Newsorted=rbind(xsorted0[1:posFC,],Separation=c("-","-","-"),xsorted0[(posFC+1):(dim(xsorted0)[1]),])
      keepOrder=c(keepOrder,posFC)
    }
    else{
      xsorted0[(posFC+1):(dim(xsorted0)[1]),] <- " "
      Newsorted=xsorted0
      keepOrder=c(keepOrder,posFC)
    }
    
    
    #output=c(output,Newsorted)
    if(nrow(output)==0){
      output=Newsorted
    }
    else{
      output=cbind(output,Newsorted)
    }
  }
  # To order clone depending on the number of sgRNA they have
  if(orderOutput){
    output=output[,(rep(order(keepOrder),each=3)*3)+rep(0:2,length(keepOrder))-2]
    
    # Since order has changed, add a space to delimitate control clone
    emptycolumn=data.frame(matrix(ncol=1,nrow=nrow(output)))
    colnames(emptycolumn)="New Clones" 
    output=cbind(output,emptycolumn)
    
  }
  
  print('Max fold change of control clones :')
  print(ALLfoldchange)
  print('The smalest fold change that will be used as threshold  :')
  print(min(ALLfoldchange))
  print('Which is control clone :')
  print(names(data)[match(min(ALLfoldchange),ALLfoldchange)+2])
  
  refFoldchange=min(ALLfoldchange)
  
  position=controlEndPos+1;
  keepOrder=c()
  outputNewClones=data.frame()
  for(i in position:length(data)){
    foldchange=c()
    x=data[,i]
    xsorted =data[order(x,decreasing = TRUE),c(1,2,i)]
    xsortedNotZscore=dataNotZscored[order(x,decreasing = TRUE),c(1,2,i)]
    if(!xsorted[1,3]==0){
      xsorted0=xsorted[xsorted[,3] != 0, ]
      
      # To find sudden drop 
      for(j in 1:maxsgRNA)
      {
        foldchange=c(foldchange,xsorted0[j,3]/xsorted0[j+1,3])
      }
      
      posFC=match(max(foldchange),foldchange)
      if(max(foldchange)>refFoldchange && xsortedNotZscore[posFC,3]>minReadCount){
        
        if(!shortOutput){
          Newsorted=rbind(xsorted0[1:posFC,],Separation=c("-","-","-"),xsorted0[(posFC+1):(dim(xsorted0)[1]),])
          keepOrder=c(keepOrder,posFC)
        }
        else{
          xsorted0[(posFC+1):(dim(xsorted0)[1]),] <- " "
          Newsorted=xsorted0
          keepOrder=c(keepOrder,posFC)
        }
        
        if(nrow(outputNewClones)==0){
          outputNewClones=Newsorted
        }
        else{
          outputNewClones=cbind(outputNewClones,Newsorted)
        }
      }
    }
  }
  
  # To order clone depending on the number of sgRNA they have
  if(orderOutput){
    outputNewClones=outputNewClones[,(rep(order(keepOrder),each=3)*3)+rep(0:2,length(keepOrder))-2]
  }
  
  output=cbind(output,outputNewClones)

  utils::write.csv(output,resultfile, row.names = FALSE,na = " ",quote=F)
}

