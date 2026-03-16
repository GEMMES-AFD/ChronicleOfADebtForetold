library(fmsb)
generate.valueMatrix<-function(ds,matrix,pars,time){
  
  # ds=data_col.ts
  # matrix=SAM
  # year=2010
  
  results=matrix(NA,nrow=nrow(matrix),ncol=ncol(matrix),dimnames = dimnames(matrix))
  period=which(ds$time==time)
  attach(ds[period,])
  attach(pars)
  for(r in 1:nrow(matrix)){
    for(c in 1:ncol(matrix)){
      cellVal<-as.character(matrix[r,c])
      if(nchar(cellVal)>0){
        #If there's an entry in the cell
        eq<-eval(parse(text=cellVal))
        results[r,c]=eq
      }
    }
  }
  detach(pars)
  detach(ds[period,])
  return(results)
}

growth<-function(var){
  return(10*(var[-1]/var[-length(var)]-1))
}


mymatplot<-function(dataset,varnames,main=NULL,location="topleft",varLegend=NULL){
  if(is.null(varLegend))
    varLegend<-varnames
  if(is.null(main))
    varLegend<-""
  dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=dataset)))
  matplot(dataset$time[1:nrow(dstemp)],dstemp, main=main ,type="l", ylab="", xlab="", lwd=2)
  if(is.numeric(varLegend)){
    if(varLegend!=-1)
      legend(location,legend=varLegend,lty=1:length(varnames),lwd=2,col=1:length(varnames),bty='n')
  }else
    legend(location,legend=varLegend,lty=1:length(varnames),lwd=2,col=1:length(varnames),bty='n')
}

mymatplotcompare<-function(datasets,varnames,location,xaxis=NA){
  dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[1]])))
  namesScen<-names(datasets)
  if(length(datasets)>1){
    for(i in 2:length(datasets)){
      dstemp=cbind(dstemp,sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[i]]))))
    }
  }
  matplot(ifelse(!is.na(xaxis),xaxis,1:nrow(dstemp)),dstemp, type="l", ylab="", xlab="", lwd=2,col=1:length(namesScen),lty=1:length(namesScen),main=varnames)
  legend(location,legend=c(namesScen),lty=1:length(namesScen),lwd=2,col=1:length(namesScen),bty='n')
}

mymatplotcompareRef<-function(datasets,varnames,location,ref=1){
  dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[1]])))
  ltys=rep(1,length(varnames))
  namesScen<-names(datasets)
  if(length(datasets)>1){
    for(i in 2:length(datasets)){
      dstemp=cbind(dstemp,sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[i]]))))
      ltys=c(ltys,rep(i,length(varnames)))
    }
  }
  dstemp=dstemp/dstemp[,ref]
  dstemp=dstemp[,-1]
  matplot(1:nrow(dstemp),dstemp, main="" ,type="l", ylab="", xlab="", lwd=2,col=1:length(varnames),lty=ltys)
  legend(location,legend=c(varnames,namesScen[-1]),lty=c(rep(1,length(varnames)),seq(1,length(namesScen))),lwd=2,col=c(1:length(varnames),rep(1,length(namesScen))),bty='n')
}

mymatplotcompareRefDif<-function(datasets,varnames,location,ref=1){
  dstemp=sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[1]])))
  ltys=rep(1,length(varnames))
  namesScen<-names(datasets)
  if(length(datasets)>1){
    for(i in 2:length(datasets)){
      dstemp=cbind(dstemp,sapply(varnames,function(x) return(eval(parse(text=x),envir=datasets[[i]]))))
      ltys=c(ltys,rep(i,length(varnames)))
    }
  }
  dstemp=dstemp-dstemp[,ref]
  dstemp=dstemp[,-1]
  matplot(1:nrow(dstemp),dstemp, main="" ,type="l", ylab="", xlab="", lwd=2,col=1:length(varnames),lty=ltys)
  legend(location,legend=c(varnames,namesScen[-1]),lty=c(rep(1,length(varnames)),seq(1,length(namesScen))),lwd=2,col=c(1:length(varnames),rep(1,length(namesScen))),bty='n')
}

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 0,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
