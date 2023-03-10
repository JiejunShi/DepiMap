library("getopt")
OptList=c('help',           'h', 0, 'logical',   'useage',
          'SampleMatrixes', 's', 1, 'character', 'Sample matrixes seperated by ",". REQUIRED.',
          'SampleNames',    'n', 1, 'character', 'Names list("," seperated) in the same length as SampleMatrixes. REQUIRED.',
          'Zmax'   ,        'z', 1, 'character', 'Zmax list("," seperated) in the length of either 1 or the same as SampleMatrixes. REQUIRED.',
          'Colors',         'c', 1, 'character', 'Color list("," seperated) in the length of either 1 or the same as SampleMatrixes. Default is blue.',
          'SmoothWindow',   'S', 1, 'integer',   'Window size(number of matrix columns) used for Smoothing. Default is 1, i.e. no smoothing.',

          'Subtrahend',     'u', 1, 'character', 'Sample Names("," seperated) which will be used as Subtrahend in Delta-Heatmaps.',
          'Minuend',        'm', 1, 'character', 'Sample Name which will be used as Minuend in Delta-Heatmaps, only 1 Minuend sample supported.',
          'DeltaZmax'   ,   'Z', 1, 'character', 'Zmax of Delta-Heatmaps("," seperated) in the length of either 1 or the same as Subtrahend Samples. REQUIRED if Subtrahend and Minuend are detected.',
          'DeltaColors',    'C', 1, 'character', '3 colors("," seperated)corresponding to -DeltaZmax,0,DeltaZmax in Delta-Heatmaps. Default is "black,white,orange".',

          'IDfiles',        'i', 1, 'character', 'ID list files("," seperated). One ID each row. Only the rows in matrixes whose id(4th column) included are presented in heatmap. Multiple horizontal panels will presented if multiple ID files detected. If absent, all rows are presented.',
          'IDorder',        'I', 0, 'logical',   'Whether sort each horizontal panel by IDfiles.',
          'WhichOrder',     'w', 1, 'character', 'Sample Names which will be used to order each horizontal panel. The length must be the same as IDfiles count. Default is the first Sample Name.',
          'OrderInterval',  'v', 1, 'character', 'Two integers indicate the columns which are used to order matrix rows. If absent, all columns are used.',
          'Decreasing',     'd', 1, 'character', 'List of 1 and/or 0 means whether order Decreasingly or not. The length must be the same as IDfiles count. Default is 1.',
          'DeltaOrder',     'o', 1, 'character', 'List of 1 and/or 0 means whether order by Delta-Heatmap or Raw heatmap. The length must be the same as IDfiles count. The Delta-Heatmap must be plotted if you want to order by it. Default is 0.',

          'LabelIDfile',    'l', 1, 'character', 'A file of IDs which will be labeled in heatmap. One ID each row.',
          'ColorResolution','r', 1, 'integer',   'Divide 1 into this number of colors. Default is 1.',
          'UseRaster',      'R', 0, 'logical',   'Whether useRaster for the heatmaps.',
          'OutTiff',        't', 1, 'character', 'Output tiff file name. Default is "Output_Heatmap.tiff".',
          'PanelWidth',     'W', 1, 'integer',   'Width of each Column Panel. Default is 300.',
          'TotalHeight',    'H', 1, 'integer',   'Total Height of the plot. Default is 1200.',
          'Xtick',          'x', 1, 'character', 'Xtick interval length indicated by either 2(TSS upstream and downstream) or 3(TSS upstream, genebody and TES downstream) integers seperated by ",". Default is 5000,5000.'
)
OptList=matrix(OptList, byrow=TRUE, ncol=5)
opt=getopt(OptList);

if(!is.null(opt$help)){
  cat(getopt(OptList, command = paste("Rscript",get_Rscript_filename()), usage=TRUE));
  q(status=1);
}
if(is.null(opt$SampleMatrixes)){
  print("[ERROR] No Input Matrix detected.");
  q(status=1)
}else{
  SampleMatrixes=unlist(strsplit(opt$SampleMatrixes,',',fixed=T))
}
if(is.null(opt$SampleNames)){
  print("[ERROR] No Sample Name detected.");
  q(status=1)
}else{
  SampleNames=unlist(strsplit(opt$SampleNames,',',fixed=T))
}
if(length(SampleMatrixes)!=length(SampleNames)){
  print("[ERROR] Count of SampleNames not equal to SampleMatrixes.");
  q(status=1)
}
if(is.null(opt$Zmax)){
  print("[ERROR] Zmax is missing.");
  q(status=1)
}else{
  Zmax=as.numeric(unlist(strsplit(opt$Zmax,',',fixed=T)))
  if(is.na(max(Zmax))){
    print("[ERROR] Zmax must be numeric.");q(status=1)
  }
  if(length(Zmax)!=1 & length(Zmax)!=length(SampleMatrixes)){
    print("[ERROR] The count of Zmax must be either 1 or the same as SampleMatrixes.");q(status=1)
  }
}
if(is.null(opt$Colors)){
  Colors="blue"
}else{
  Colors=unlist(strsplit(opt$Colors,',',fixed=T))
  if(length(Colors)!=1 & length(Colors)!=length(SampleMatrixes)){
    print("[ERROR] The count of Colors must be either 1 or the same as SampleMatrixes.");q(status=1)
  }
}
if(is.null(opt$SmoothWindow)){
  SmoothWindow=1
}else{
  SmoothWindow=as.integer(opt$SmoothWindow)
}
if(is.null(opt$Subtrahend)){
  Subtrahend=NA
}else{
  Subtrahend=unlist(strsplit(opt$Subtrahend,',',fixed=T))
  if(length(setdiff(Subtrahend,SampleNames))>0){
    print("[ERROR] All Subtrahend samples must be in SampleNames.");q(status=1)
  }
}
if(is.null(opt$Minuend)){
  Minuend=NA
}else{
  Minuend=as.character(opt$Minuend)
  if(!Minuend%in%SampleNames){
    print("[ERROR] The Minuend sample must be in SampleNames.");q(status=1)
  }
}
if(is.na(Subtrahend) & (!is.na(Minuend))){
  print("[ERROR] The Minuend is detected but Subtrahend is missing.");
  q(status=1)
}
if(is.na(Minuend) & (!is.na(Subtrahend))){
  print("[ERROR] The Subtrahend is detected but Minuend is missing.");
  q(status=1)
}
if((!is.na(Subtrahend)) & (!is.na(Minuend))){
  if(is.null(opt$DeltaZmax)){
    print("[ERROR] DeltaZmax is missing.");q(status=1)
  }else{
    DeltaZmax=as.numeric(unlist(strsplit(opt$DeltaZmax,',',fixed=T)))
    if(is.na(max(DeltaZmax))){
      print("[ERROR] DeltaZmax must be numeric.");q(status=1)
    }
    if(length(DeltaZmax)!=1 & length(DeltaZmax)!=length(Subtrahend)){
      print("[ERROR] The count of DeltaZmax must be either 1 or the same as Subtrahend samples.");q(status=1)
    }
  };
  if(is.null(opt$DeltaColors)){
    DeltaColors=c("black","white","orange")
  }else{
    DeltaColors=unlist(strsplit(opt$DeltaColors,',',fixed=T))
    if(length(DeltaColors)!=3){
      print("[ERROR] Count of DeltaColors must be 3.");q(status=1)
    }
  }
}
if(is.null(opt$IDfiles)){
  IDfiles=NA
}else{
  IDfiles=unlist(strsplit(opt$IDfiles,',',fixed=T))
}
if(is.null(opt$IDorder)){
  IDorder=FALSE
}else{
  IDorder=TRUE
}
if(is.null(opt$WhichOrder)){
  WhichOrder=SampleNames[1]
}else{
  WhichOrder=unlist(strsplit(opt$WhichOrder,',',fixed=T))
  if(length(WhichOrder)!=length(IDfiles)){
    print("[ERROR] The length of WhichOrder must be the same as IDfiles count.");q(status=1)
  }
  if(length(setdiff(WhichOrder,SampleNames))>0){
    print("[ERROR] The WhichOrder must be in SampleNames.");q(status=1)
  }
}
if(is.null(opt$OrderInterval)){
  OrderInterval=NA
}else{
  OrderInterval=as.integer(unlist(strsplit(opt$OrderInterval,',',fixed=T)))
  if(is.na(max(OrderInterval))|length(OrderInterval)!=2){
    print("[ERROR] If provided, OrderInterval must be 2 integers.");q(status=1)
  }
  OrderInterval=sort(OrderInterval)
}
if(is.null(opt$Decreasing)){
  Decreasing=1
}else{
  Decreasing=as.numeric(unlist(strsplit(opt$Decreasing,',',fixed=T)))
  if(length(setdiff(Decreasing,c(0,1)))>0){
    print("[ERROR] The Decreasing must be a list of 0 and/or 1.");q(status=1)
  }
  if(length(Decreasing)!=length(IDfiles)){
    print("[ERROR] The length of Decreasing must be the same as IDfiles count.");q(status=1)
  }
}
if(is.null(opt$DeltaOrder)){
  DeltaOrder=0
}else{
  DeltaOrder=as.numeric(unlist(strsplit(opt$DeltaOrder,',',fixed=T)))
  if(length(setdiff(DeltaOrder,c(0,1)))>0){
    print("[ERROR] The DeltaOrder must be a list of 0 and/or 1.");q(status=1)
  }
  if(length(DeltaOrder)!=length(IDfiles)){
    print("[ERROR] The length of DeltaOrder must be the same as IDfiles count.");q(status=1)
  }
  if(length(setdiff(WhichOrder[DeltaOrder==1],Subtrahend))>0){
    print('[ERROR] The Delta heatmap must be plotted if you want to order by it. Please check whether the samples with DeltaOrder=1 are all in "Subtrahend".');q(status=1)
  }
}
if(is.null(opt$LabelIDfile)){
  LabelID=NA
}else{
  LabelID=as.character(opt$LabelIDfile)
}
if(is.null(opt$ColorResolution)){
  ColorResolution=1
}else{
  ColorResolution=as.integer(opt$ColorResolution)
}
if(is.null(opt$UseRaster)){
  UseRaster=FALSE
}else{
  UseRaster=TRUE
}
if(is.null(opt$OutTiff)){
  OutTiff="Output_Heatmap.tiff"
}else{
  OutTiff=as.character(opt$OutTiff)
}
if(is.null(opt$PanelWidth)){
  PanelWidth=300
}else{
  PanelWidth=as.numeric(opt$PanelWidth)
}
if(is.null(opt$TotalHeight)){
  TotalHeight=1200
}else{
  TotalHeight=as.numeric(opt$TotalHeight)
}
if(is.null(opt$Xtick)){
  Xtick=c(5000,5000)
}else{
  Xtick=as.integer(unlist(strsplit(opt$Xtick,',',fixed=T)))
  if(is.na(max(Xtick))|(length(Xtick)!=2 & length(Xtick)!=3)){
    print("[ERROR] Xtick must be either 2 or 3 integers.");q(status=1)
  }
}
##################################
# Functions read/order matrix
##################################
ReadMatrix<-function(MatrixFileName,IDlist){
  aMatrix=read.delim(gzfile(MatrixFileName),header=F,skip=1)
  rownames(aMatrix)=as.character(aMatrix$V4);
  aMatrix=aMatrix[,7:ncol(aMatrix)];
  if(!is.na(IDlist)){
    IDlist=as.character(read.delim(IDlist,header = F)$V1)
    aMatrix=aMatrix[IDlist,]
  }
  aMatrix=aMatrix[order(rownames(aMatrix)),]
  return(aMatrix)
}
HorizontalPanelOrder<-function(OrderSample,WhetherDecreasing,WhetherDeltaOrder,IDlist,WhetherIDorder){
  if(WhetherIDorder==TRUE){
    aOrder=rev(as.character(read.delim(IDlist,header = F)$V1))
  }else{
    if(WhetherDecreasing==0){
      decreasing_flag=TRUE
    }else{
      decreasing_flag=FALSE
    }
    mat1=ReadMatrix(MatrixFileName=SampleMatrixes[SampleNames==OrderSample],IDlist)
    if(!is.na(OrderInterval[1])){
      OrderColumns=sort(intersect(1:ncol(mat1),OrderInterval[1]:OrderInterval[2]))
      if(length(OrderColumns)==0){
        OrderColumns=NA
        print(paste("Both OrderInterval exceed matrix column bondary(1:",ncol(mat1),")! Use all columns to order as default.",sep=""))
      }else{
        print(paste("Use ",OrderColumns[1],":",OrderColumns[length(OrderColumns)]," columns to order.",sep=""))
      }
    }else{
      OrderColumns=NA
    }
    if(WhetherDeltaOrder==1){
      mat2=ReadMatrix(MatrixFileName=SampleMatrixes[SampleNames==Minuend],IDlist)
      mat_delta=mat1-mat2
      if(!is.na(OrderColumns[1])){mat_delta=mat_delta[,OrderColumns]}
      mat_delta=mat_delta[order(rowSums(mat_delta),decreasing=decreasing_flag),]
      aOrder=rownames(mat_delta)
    }else{
      if(!is.na(OrderColumns[1])){mat1=mat1[,OrderColumns]}
      mat1=mat1[order(rowSums(mat1),decreasing=decreasing_flag),]
      aOrder=rownames(mat1)
    }
  };
  return(aOrder)
}
AdjacentAverage<-function(x=c(),w=1){
  x1=c();
  for(i in 1:length(x)){
    x0=mean(x[max(1,(i-round(w/2))):min((i+round(w/2)),length(x))]);
    x1=c(x1,x0);
  };
  return(x1)
}
##################################
# Pre-process
##################################
if(is.na(Subtrahend)){
  VerticalPanel=length(SampleNames)+2
}else{
  VerticalPanel=length(SampleNames)+length(Subtrahend)+2
}
if(length(Decreasing)==1 & length(IDfiles)>1){
  Decreasing=rep(Decreasing,length(IDfiles))
}
if(length(DeltaOrder)==1 & length(IDfiles)>1){
  DeltaOrder=rep(DeltaOrder,length(IDfiles))
}
TotalRowCount=c()
for(i in 1:length(IDfiles)){
  aRowOrder=HorizontalPanelOrder(OrderSample=WhichOrder[i],WhetherDecreasing=Decreasing[i],WhetherDeltaOrder=DeltaOrder[i],IDlist=IDfiles[i],WhetherIDorder=IDorder)
  assign(paste("Panel",i,"RowOrder",sep=""),aRowOrder)
  TotalRowCount=c(TotalRowCount,length(aRowOrder))
}
Bottom=60/TotalHeight;Top=1-60/TotalHeight;
##################################
# Plot
##################################
tiff(OutTiff,width=PanelWidth*VerticalPanel,height=TotalHeight,compression="lzw",res=600)
for(s in 1:length(SampleNames)){
  if(length(Colors)==1){
    aColor=Colors
  }else{
    aColor=Colors[s]
  }
  if(length(Zmax)==1){
    aZmax=Zmax
  }else{
    aZmax=Zmax[s]
  }
  aRawColor=colorRampPalette(c("white",aColor))(aZmax*ColorResolution+1)
  x1=s/VerticalPanel;
  x2=(s+1)/VerticalPanel;
  for(i in 1:length(IDfiles)){
    aMatrix=ReadMatrix(MatrixFileName=SampleMatrixes[s],IDlist=IDfiles[i])
    aMatrix=aMatrix[get(paste("Panel",i,"RowOrder",sep="")),]
    if(SmoothWindow>1){
      aMatrix=apply(aMatrix,1,AdjacentAverage,w=SmoothWindow);aMatrix=t(aMatrix);
    }
    aMatrix[which(aMatrix>aZmax,arr.ind = T)]=aZmax

    if(i==1){
      y2=Top
    }else{
      y2=Top-sum(TotalRowCount[1:(i-1)])/sum(TotalRowCount)*(Top-Bottom)
    }
    y1=y2-TotalRowCount[i]/sum(TotalRowCount)*(Top-Bottom);
    par(fig=c(x1,x2,y1,y2),mai=c(0,0,10/TotalHeight,10/TotalHeight),new=T,mgp=c(1,0.1,0))
    image(z=t(aMatrix),col=aRawColor[0:ceiling(max(aMatrix,na.rm=T)*ColorResolution)+1],useRaster=UseRaster,xaxt="n",yaxt="n",bty="n");box(lwd=0.2)
    if(s==1){axis(2,at=0.5,labels=paste("n=",nrow(aMatrix),sep=""),tick=F,cex.axis=0.15,line=0)}
    if(i==length(IDfiles)){
      if(length(Xtick)==2){
        axis(1,at=c(0,Xtick[1]/sum(Xtick),1),labels=c(paste(-Xtick[1]/1000,"kb",sep=""),0,paste(Xtick[2]/1000,"kb",sep="")),las=2,tick=T,cex.axis=0.15,line=0,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
      }else if(length(Xtick)==3){
        axis(1,at=c(0,Xtick[1]/sum(Xtick),sum(Xtick[1:2])/sum(Xtick),1),labels=c(paste(-Xtick[1]/1000,"kb",sep=""),"TSS","TES",paste(Xtick[3]/1000,"kb",sep="")),las=2,tick=T,cex.axis=0.15,line=0,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
      }
    }
    if((!is.na(LabelID)) & s==(VerticalPanel-2)){
      ToLabel=intersect(as.character(read.delim(LabelID,header = F)$V1),rownames(aMatrix));
      ToLabel=which(rownames(aMatrix)%in%ToLabel)
      axis(4,at=(ToLabel-1)/(nrow(aMatrix)-1),labels=NA,las=2,tick=T,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
      axis(4,at=(1:length(ToLabel))/length(ToLabel),labels=rownames(aMatrix)[ToLabel],las=2,tick=F,cex.axis=0.15,line=0.25)
      segments(x0=rep(1.05,length(ToLabel)),y0=(ToLabel-1)/(nrow(aMatrix)-1),x1=rep(1.2,length(ToLabel)),y1=(1:length(ToLabel))/length(ToLabel),lty=1,lwd=0.2,xpd=NA)
    }
  }
  par(fig=c(x1,x2,Top+(1-Top)/4,Top+(1-Top)/2),mai=c(0,0,0,10/TotalHeight),new=T,mgp=c(1,0.1,0))
  image(z=matrix(0:(aZmax*ColorResolution),ncol=1),col=aRawColor,xaxt="n",yaxt="n",bty="n");box(lwd=0.2)
  axis(1,at=0,labels=0,tick=F,cex.axis=0.15,hadj=0,line=-0.75);
  axis(1,at=1,labels=aZmax,tick=F,cex.axis=0.15,hadj=1,line=-0.75)
  axis(3,at=0.5,labels=SampleNames[s],tick=F,cex.axis=0.15,line=-0.2)
}
if(!is.na(Subtrahend)){
  for(s in 1:length(Subtrahend)){
    if(length(DeltaZmax)==1){
      aDeltaZmax=DeltaZmax
    }else{
      aDeltaZmax=DeltaZmax[s]
    }
    aDeltaColor=colorRampPalette(DeltaColors)(2*aDeltaZmax*ColorResolution+1)
    x1=(s+length(SampleNames))/VerticalPanel;
    x2=(s+length(SampleNames)+1)/VerticalPanel;
    for(i in 1:length(IDfiles)){
      mat1=ReadMatrix(MatrixFileName=SampleMatrixes[SampleNames==Subtrahend[s]],IDlist=IDfiles[i])
      mat1=mat1[get(paste("Panel",i,"RowOrder",sep="")),]
      mat2=ReadMatrix(MatrixFileName=SampleMatrixes[SampleNames==Minuend],IDlist=IDfiles[i])
      mat2=mat2[get(paste("Panel",i,"RowOrder",sep="")),]
      aMatrix=mat1-mat2
      if(SmoothWindow>1){
        aMatrix=apply(aMatrix,1,AdjacentAverage,w=SmoothWindow);aMatrix=t(aMatrix);
      }
      aMatrix[which(aMatrix>aDeltaZmax,arr.ind = T)]=aDeltaZmax
      aMatrix[which(aMatrix<(-aDeltaZmax),arr.ind = T)]=(-aDeltaZmax)

      if(i==1){
        y2=Top
      }else{
        y2=Top-sum(TotalRowCount[1:(i-1)])/sum(TotalRowCount)*(Top-Bottom)
      }
      y1=y2-TotalRowCount[i]/sum(TotalRowCount)*(Top-Bottom);
      par(fig=c(x1,x2,y1,y2),mai=c(0,0,10/TotalHeight,10/TotalHeight),new=T,mgp=c(1,0.1,0))
      image(z=t(aMatrix),col=aDeltaColor[floor(min(aMatrix,na.rm=T)*ColorResolution):ceiling(max(aMatrix,na.rm=T)*ColorResolution)+aDeltaZmax*ColorResolution+1],useRaster=UseRaster,xaxt="n",yaxt="n",bty="n");box(lwd=0.2)
      if(i==length(IDfiles)){
        if(length(Xtick)==2){
          axis(1,at=c(0,Xtick[1]/sum(Xtick),1),labels=c(paste(-Xtick[1]/1000,"kb",sep=""),0,paste(Xtick[2]/1000,"kb",sep="")),las=2,tick=T,cex.axis=0.15,line=0,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
        }else if(length(Xtick)==3){
          axis(1,at=c(0,Xtick[1]/sum(Xtick),sum(Xtick[1:2])/sum(Xtick),1),labels=c(paste(-Xtick[1]/1000,"kb",sep=""),"TSS","TES",paste(Xtick[3]/1000,"kb",sep="")),las=2,tick=T,cex.axis=0.15,line=0,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
        }
      }
      if((!is.na(LabelID)) & s==(VerticalPanel-length(SampleNames)-2)){
        ToLabel=intersect(as.character(read.delim(LabelID,header = F)$V1),rownames(aMatrix));
        ToLabel=which(rownames(aMatrix)%in%ToLabel)
        axis(4,at=(ToLabel-1)/(nrow(aMatrix)-1),labels=NA,las=2,tick=T,lwd=0.1,lwd.ticks=0.2,tck=-0.05)
        axis(4,at=(1:length(ToLabel))/length(ToLabel),labels=rownames(aMatrix)[ToLabel],las=2,tick=F,cex.axis=0.15,line=0.25)
        segments(x0=rep(1.05,length(ToLabel)),y0=(ToLabel-1)/(nrow(aMatrix)-1),x1=rep(1.2,length(ToLabel)),y1=(1:length(ToLabel))/length(ToLabel),lty=1,lwd=0.2,xpd=NA)
      }
    }
    par(fig=c(x1,x2,Top+(1-Top)/4,Top+(1-Top)/2),mai=c(0,0,0,10/TotalHeight),new=T,mgp=c(1,0.1,0))
    image(z=matrix(0:(aDeltaZmax*ColorResolution*2),ncol=1),col=aDeltaColor,xaxt="n",yaxt="n",bty="n");box(lwd=0.2)
    axis(1,at=0,labels=(-aDeltaZmax),tick=F,cex.axis=0.15,hadj=0,line=-0.75);
    axis(1,at=0.5,labels=0,tick=F,cex.axis=0.15,hadj=0.5,line=-0.75);
    axis(1,at=1,labels=aDeltaZmax,tick=F,cex.axis=0.15,hadj=1,line=-0.75)
    axis(3,at=0.5,labels=paste(Subtrahend[s],Minuend,sep="-"),tick=F,cex.axis=0.15,line=-0.2)
  }
}
dev.off()
