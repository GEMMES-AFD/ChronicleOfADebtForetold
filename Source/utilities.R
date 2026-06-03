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
ggparcoord_ind_yaxis <- function(
    data,
    truth = NULL, 
    truthPointSize = 2, 
    columns = 1:ncol(data),
    groupColumn = NULL, 
    alphaLines = 1, 
    nbreaks = 4, 
    axis_font_size = 3, 
    linewidth = 1,
    custom_breaks = NULL,
    axis_normal = "none",
    axis_normal_coeff = 1
) {
  
  # select the variables to plot
  data_subset <- data %>% dplyr::select(columns)
  if (!is.null(custom_breaks)){
    custom_breaks_subset <- custom_breaks %>% dplyr::select(columns)
  }
  # re-order truth to match columns
  col_names <- data_subset %>% names
  if (!is.null(truth)) {
    truth <- truth %>% select(col_names)
    data_subset <- data_subset %>% rbind(truth)
  } 
  
  # Calculate the axis breaks for each variable on the *original* scale.
  # Note that the breaks computed by pretty() are guaranteed to contain all of 
  # the data. We include truth in these breaks, just in case one of the true 
  # points falls outside the range of the data (can easily happen in the context
  # of comparing parameter estimates to the true values).
  
  if (is.null(custom_breaks)){
    breaks_df <- data_subset %>% 
      stack %>%           # convert to long format
      group_by(ind) %>%   # group by the plotting variables
      summarize(breaks = pretty(values, n = nbreaks))
  } else  {
    breaks_df <- custom_breaks_subset %>% 
      stack %>%           # convert to long format
      group_by(ind) %>%   # group by the plotting variables
      summarize(breaks = pretty(values, n = nbreaks)) 
  }
  
  # Normalise the breaks to be between 0 and 1, and set the coordinates of the 
  # tick marks. Importantly, if we want the axis heights to be the same, the 
  # breaks need to be normalised to be between exactly 0 and 1. 
  axis_df <- breaks_df %>% 
    mutate(yval = (breaks - min(breaks))/(max(breaks) - min(breaks))) %>%
    mutate(xmin = as.numeric(ind) - 0.05, 
           xmax = as.numeric(ind),
           x_text = as.numeric(ind) - 0.2)
  
  
  axis_df <- axis_df %>% mutate(breaks2 = ifelse(ind %in% axis_normal, breaks*axis_normal_coeff, breaks))
  
  # Calculate the co-ordinates for our axis lines:
  axis_line_df <- axis_df %>% 
    group_by(ind) %>%
    summarize(min = min(yval), max = max(yval))
  
  # Getting the minimum/maximum breaks on the original scale, to scale the 
  # data in the same manner that we scaled the breaks
  minmax_breaks <- breaks_df %>%
    summarize(min_break = min(breaks), max_break = max(breaks)) %>% 
    tibble::column_to_rownames(var = "ind")
  
  # Normalise the original data in the same way that the breaks were normalised.
  # This ensures that the scaling is correct. 
  # Do the same for the truth points, if they exist.
  lines_df <- data %>% select(columns) 
  for (col in col_names) {
    lines_df[, col] <- (lines_df[, col] - minmax_breaks[col, "min_break"]) / ( minmax_breaks[col, "max_break"] -  minmax_breaks[col, "min_break"])
    if (!is.null(truth)) {
      truth[, col] <- (truth[, col] - minmax_breaks[col, "min_break"]) / ( minmax_breaks[col, "max_break"] -  minmax_breaks[col, "min_break"])
    }
  }
  
  # Reshape original data (and truth):
  lines_df <- lines_df %>%
    mutate(row = row_number()) %>% # need row information to group individual rows
    bind_cols(data[, groupColumn, drop = FALSE]) %>% # need groupColumn for colour aesthetic
    reshape2::melt(id.vars = c("row", groupColumn), 
                   # choose names that are consistent with stack() above:
                   value.name = "values", variable.name = "ind") 
  
  # Reshape truth, as above
  if (!is.null(truth)) {
    truth <- truth %>%
      mutate(row = row_number()) %>% # need row information to group individual rows
      reshape2::melt(id.vars = c("row"), 
                     # choose names that are consistent with stack():
                     value.name = "values", variable.name = "ind") 
  }
  
  if (!is.numeric(linewidth)){
    # Now plot: 
    gg <- ggplot() + 
      geom_line(data = lines_df %>% sample_n(nrow(.)), # permute rows to prevent one group dominating over another
                aes_string(x = "ind", y = "values", group = "row", colour = groupColumn, alpha=alphaLines, linewidth = linewidth)) +
      geom_segment(data = axis_line_df, aes(x = ind, xend = ind, y = min, yend = max),
                   inherit.aes = FALSE) +
      geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = yval, yend = yval),
                   inherit.aes = FALSE) +
      geom_text(data = axis_df, aes(x = x_text, y = yval, label = breaks2),
                inherit.aes = FALSE, size = axis_font_size) 
    
    if (!is.null(truth)) {
      gg <- gg + geom_point(data = truth, aes(x = ind, y = values), 
                            inherit.aes = FALSE, colour = "red", size = truthPointSize)
    }
  } else {
    gg <- ggplot() + 
      geom_line(data = lines_df %>% sample_n(nrow(.)), # permute rows to prevent one group dominating over another
                aes_string(x = "ind", y = "values", group = "row", colour = groupColumn, alpha=alphaLines), linewidth = linewidth)  +
      geom_segment(data = axis_line_df, aes(x = ind, xend = ind, y = min, yend = max),
                   inherit.aes = FALSE) +
      geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = yval, yend = yval),
                   inherit.aes = FALSE) +
      geom_text(data = axis_df, aes(x = x_text, y = yval, label = breaks),
                inherit.aes = FALSE, size = axis_font_size) 
    
    if (!is.null(truth)) {
      gg <- gg + geom_point(data = truth, aes(x = ind, y = values), 
                            inherit.aes = FALSE, colour = "red", size = truthPointSize) 
    }
  }
  gg <- gg + theme_bw() + 
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(),
          axis.ticks =  element_blank(), 
          axis.text.y = element_blank()) 
  
  
  
  return(gg)
}



ggparcoord_ind_yaxis2 <- function(
    data,
    truth = NULL, 
    truthPointSize = 2, 
    columns = 1:ncol(data),
    groupColumn = NULL, 
    alphaLines = 1, 
    nbreaks = 4, 
    axis_font_size = 3,
    main =NULL,
    linewidth = 1,
    varlabs = NULL,
    palette = "viridis",
    direction =1,
    begin = 0,
    end = 1,
    sd_var = NULL,
    min_width = 0.1,
    max_width = 0.5,
    width_label = "Standard \n deviation \n of outcome",
    threshold = 1,
    alpha = F,
    min_alpha = 0.5,
    max_alpha = 1,
    color_label = "Mean \n value \n of outcome",
    axis_width = 1,
    guideinstruction_mean ="legend",
    guideinstruction_var = "legend",
    size_xticks = 10,
    custom_breaks = NULL,
    round = 0,
    custom_lower_break_outcome = NULL,
    custom_higher_break_outcome = NULL,
    logoutcome=F,
    highlight = F,
    minhighlight = 0.9,
    maxhighlight = 1,
    col_par=2
) {
  
  require(gghighlight)
  # select the variables to plot
  if (logoutcome){
    data <- data%>%
      mutate_at(groupColumn, logzero)
    data <- data %>%
      mutate_at(sd_var, logzero)
  }
  
  data_subset <- data %>% dplyr::select(columns)
  if (!is.null(custom_breaks)){
    if (logoutcome){
      custom_breaks <- custom_breaks%>%
        mutate_at(groupColumn, logzero)
      custom_breaks <- custom_breaks %>%
        mutate_at(sd_var, logzero)
    }
    custom_breaks_subset <- custom_breaks %>% dplyr::select(columns)
  }
  
  
  # re-order truth to match columns
  col_names <- data_subset %>% names
  if (!is.null(truth)) {
    truth <- truth %>% dplyr::select(col_names)
    data_subset <- data_subset %>% rbind(truth)
  } 
  
  # Calculate the axis breaks for each variable on the *original* scale.
  # Note that the breaks computed by pretty() are guaranteed to contain all of 
  # the data. We include truth in these breaks, just in case one of the true 
  # points falls outside the range of the data (can easily happen in the context
  # of comparing parameter estimates to the true values).
  
  
  
  if (is.null(custom_breaks)){
    breaks_df <- data_subset %>% 
      stack %>%           # convert to long format
      group_by(ind) %>%   # group by the plotting variables
      summarize(breaks = pretty(values, n = nbreaks))
  } else  {
    breaks_df <- custom_breaks_subset %>% 
      stack %>%           # convert to long format
      group_by(ind) %>%   # group by the plotting variables
      summarize(breaks = pretty(values, n = nbreaks)) 
  }
  
  if (!is.null(custom_lower_break_outcome)){
    index <- which(breaks_df$ind == groupColumn)
    index_min <- min(index) -1 + which(breaks_df[index,2] == min(breaks_df[index,2]))
    breaks_df$breaks[index_min] <- custom_lower_break_outcome
    breaks_df <- breaks_df[-(index_min+1),]
  }
  
  if (!is.null(custom_higher_break_outcome)){
    index <- which(breaks_df$ind == groupColumn)
    index_max <- min(index) -1 + which(breaks_df[index,2] == max(breaks_df[index,2]))
    breaks_df$breaks[index_max] <- custom_higher_break_outcome
    breaks_df <- breaks_df[-(index_max-1),]
  }
  # Normalise the breaks to be between 0 and 1, and set the coordinates of the 
  # tick marks. Importantly, if we want the axis heights to be the same, the 
  # breaks need to be normalised to be between exactly 0 and 1. 
  axis_df <- breaks_df %>% 
    mutate(yval = (breaks - min(breaks))/(max(breaks) - min(breaks))) %>%
    mutate(xmin = as.numeric(ind) - 0.05, 
           xmax = as.numeric(ind),
           x_text = as.numeric(ind) - 0.2)
  
  
  
  # Calculate the co-ordinates for our axis lines:
  axis_line_df <- axis_df %>% 
    group_by(ind) %>%
    summarize(min = min(yval), max = max(yval))
  
  # Getting the minimum/maximum breaks on the original scale, to scale the 
  # data in the same manner that we scaled the breaks
  minmax_breaks <- breaks_df %>%
    summarize(min_break = min(breaks), max_break = max(breaks)) %>% 
    tibble::column_to_rownames(var = "ind")
  
  qmeanvar <-  (minmax_breaks[groupColumn, "max_break"] -  minmax_breaks[groupColumn, "min_break"])*seq(0,1,0.2) + minmax_breaks[groupColumn, "min_break"]
  #threshold <-  (minmax_breaks[groupColumn, "max_break"] -  minmax_breaks[groupColumn, "min_break"])*threshold + minmax_breaks[groupColumn, "min_break"]
  
  # Normalise the original data in the same way that the breaks were normalised.
  # This ensures that the scaling is correct. 
  # Do the same for the truth points, if they exist.
  lines_df <- data %>% dplyr::select(columns) 
  for (col in col_names) {
    lines_df[, col] <- (lines_df[, col] - minmax_breaks[col, "min_break"])/(minmax_breaks[col, "max_break"] -  minmax_breaks[col, "min_break"])
    if (!is.null(truth)) {
      truth[, col] <- (truth[, col] - minmax_breaks[col, "min_break"])/(minmax_breaks[col, "max_break"] -  minmax_breaks[col, "min_break"])
    }
  }
  
  # Reshape original data (and truth):
  lines_df <- lines_df %>%
    mutate(row = row_number()) %>% # need row information to group individual rows
    reshape2::melt(id.vars = c("row"), 
                   # choose names that are consistent with stack() above:
                   value.name = "values", variable.name = "ind") 
  
  # Reshape truth, as above
  if (!is.null(truth)) {
    truth <- truth %>%
      mutate(row = row_number()) %>% # need row information to group individual rows
      reshape2::melt(id.vars = c("row"), 
                     # choose names that are consistent with stack():
                     value.name = "values", variable.name = "ind") 
  }
  
  if (!is.null(groupColumn)){
    lines_df_groupcolumn <- lines_df %>%
      filter(ind == groupColumn)
    lines_df <- lines_df %>%
      mutate(groupColumn2 = NA)
    for (i in 1:nrow(lines_df_groupcolumn)){
      lines_df$groupColumn2[i + seq(0,(ncol(data_subset)-1)*(nrow(lines_df_groupcolumn)), max(lines_df$row))] <- lines_df_groupcolumn$values[i]
    }
    if (direction < 0){
      lines_df <- lines_df %>%
        mutate(groupColumn3 = ifelse(groupColumn2 < threshold, 0, groupColumn2), 
               groupColumn4 = (groupColumn2 - min(groupColumn2))/(max(groupColumn2)-min(groupColumn2)))
    } else {
      lines_df <- lines_df %>%
        mutate(groupColumn3 = ifelse(groupColumn2 < threshold, 0, groupColumn2),
               groupColumn4 = (groupColumn2 - min(groupColumn2))/(max(groupColumn2)-min(groupColumn2)))
    }
  }
  
  
  if (!is.null(varlabs)){
    if (length(varlabs) == length(unique(axis_df$ind))){
      levels(lines_df$ind) <- varlabs
      levels(axis_line_df$ind) <- varlabs
    } else {
      print("Error: Not enough varlabs")
    }
  }
  
  if (is.null(sd_var)){
    lines_df <- lines_df %>%
      mutate(sd_var = NA)
    lines_df$sd_var <- 1
  } else {
    sd_var2 <- min_width + max_width*(data[, sd_var] - min(data[, sd_var]))/(max(data[, sd_var]) - min(data[, sd_var]))
    lines_df <- lines_df %>%
      mutate(sd_var = NA)
    for (i in 1:length(sd_var2[,1])){
      lines_df$sd_var[i + seq(0,(ncol(data_subset)-1)*(length(sd_var2[,1])), max(lines_df$row))] <- sd_var2[i,1] 
    }
  }
  
  
  if (highlight){
    minhighlight = quantile(lines_df$groupColumn2, minhighlight, na.rm=T)
    maxhighlight = quantile(lines_df$groupColumn2, maxhighlight, na.rm=T)
  } else {
    minhighlight = 0.99*min(lines_df$groupColumn2)
    maxhighlight = 1.01*max(lines_df$groupColumn2)
  }
  if (alpha){
    alpha <- 'groupColumn2'
  } else {
    alpha <- 1
  }
  
  if (length(unique(lines_df$sd_var)) >1){
    qsdvar <- round((max(custom_breaks[,sd_var]) -  min(custom_breaks[,sd_var]))*seq(0,1,0.2) + min(custom_breaks[,sd_var]),round)
    
    gg <- ggplot(data = lines_df, 
                 aes_string(x = "ind", y = "values", group = "row", colour = paste("groupColumn", col_par, sep=""))) +
      ggtitle(main)+
      geom_borderline(stat='smooth', aes_string(alpha = alpha), method = "loess", bordercolor = "white") +
      gghighlight(groupColumn2 > minhighlight & groupColumn2 <= maxhighlight)+
      scale_linewidth(width_label, range=c(0.1,2), breaks = seq(0,1,0.2), labels = qsdvar, guide = "none")+
      scale_color_viridis(color_label, option = palette, direction = direction, begin = begin, end = end, guide = guideinstruction_mean, limits=c(0,1), breaks = seq(0,1,0.2), labels = qmeanvar)+
      scale_alpha(range=c(min_alpha,max_alpha), guide = "none")+
      geom_segment(data = axis_line_df, aes(x = ind, xend = ind, y = min, yend = max), linewidth = axis_width,
                   inherit.aes = FALSE) +
      geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = yval, yend = yval),
                   linewidth = axis_width,
                   inherit.aes = FALSE) +
      geom_text(data = axis_df, aes(x = x_text, y = yval, label = breaks),
                inherit.aes = FALSE, size = axis_font_size)
    
  } else {
    
    gg <- ggplot(data = lines_df, 
                 aes_string(x = "ind", y = "values", group = "row", colour = paste("groupColumn", col_par, sep=""))) +
      ggtitle(main)+
      geom_borderline(stat='smooth',  aes_string(alpha = alpha),  method = "loess", linewidth =linewidth, bordercolor = "white") +
      gghighlight(groupColumn2 > minhighlight & groupColumn2 <= maxhighlight)+
      scale_linewidth(width_label, range = c(0.1,1), guide = guideinstruction_var)+
      scale_color_viridis(color_label, option = palette, direction = direction, begin = begin, end = end, guide = "none", limits=c(0,1), breaks = seq(0,1,0.2), labels = qmeanvar)+
      scale_alpha(range=c(min_alpha,max_alpha), guide = "none")+
      geom_segment(data = axis_line_df, aes(x = ind, xend = ind, y = min, yend = max), linewidth = axis_width,
                   inherit.aes = FALSE) +
      geom_segment(data = axis_df, aes(x = xmin, xend = xmax, y = yval, yend = yval),
                   linewidth = axis_width,
                   inherit.aes = FALSE) +
      geom_text(data = axis_df, aes(x = x_text, y = yval, label = breaks),
                inherit.aes = FALSE, size = axis_font_size, fontface = "bold")
  }
  
  
  if (!is.null(truth)) {
    gg <- gg + geom_point(data = truth, aes(x = ind, y = values), 
                          inherit.aes = FALSE, colour = "red", size = truthPointSize)
  }
  
  gg <- gg + theme_bw() + 
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(),
          axis.ticks =  element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = size_xticks))
  
  return(gg)
}

