plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", color = column)) +
    geom_point() +
    facet_wrap(~year) +
    coord_fixed() + theme_bw()
}

# https://stackoverflow.com/questions/8261590/write-list-to-a-text-file-preserving-names-r
save_sanity_text <- function(x, fil){ 
  nams=names(x) 
  cat(nams[1], ": ",  x[[1]], "\n", file=fil, append=FALSE)
  for (i in 2:length(x) ){ cat(nams[i], ": ",  x[[i]], "\n", file=fil, append=TRUE) }
}

calc_AUC <- function(pred, obs){
  obs <- as.numeric(as.character(obs))
  predict = ROCR::prediction(as.vector(pred), obs)
  AUC <- round(unlist(slot(performance(predict,"auc"),"y.values")),3)
  return(AUC)
}
calc_RMSE <- function(pred, obs){
  RMSE <- round(sqrt(mean((pred-obs)^2)),3)
  return(RMSE)
}

# total (site) = 6*fyke + 2*trammel45 + 2*trammel37
# sum over sites to get index by year and age
get_total_pred_byyear <- function(df, pred){
  df$pred <- pred
  df[df$gear == 'Trammel (45 mm, garn)','pred'] = 2 * df[df$gear == 'Trammel (45 mm, garn)','pred']
  df[df$gear == 'Trammel (37 mm, garn)','pred'] = 2 * df[df$gear == 'Trammel (37 mm, garn)','pred']
  df[df$gear == 'Fyke (ruse)','pred'] = 6 * df[df$gear == 'Fyke (ruse)','pred']
  bysite <- df %>% group_by(year,site) %>% summarize(total=sum(pred)) %>% as.data.frame
  ind <- bysite %>% group_by(year) %>% summarize(total=sum(total)) %>% pull(total)
  return(ind)
}

check_boot <- function(df, type='set'){
  if(type=='set'){
    CAA <- df %>% group_by(year) %>% summarize(catch_0=sum(catch_0),
                                           catch_1=sum(catch_1),
                                           catch_2=sum(catch_2),
                                           catch_3=sum(catch_3),
                                           catch_4=sum(catch_4),
                                           catch_5=sum(catch_5),
                                           catch_6=sum(catch_6),
                                           catch_7=sum(catch_7)) %>% select(-year)
  }
  if(type == 'site'){
    CAA <- df %>% group_by(year, Age) %>% summarize(catch=sum(Catch)) %>% 
                  pivot_wider(names_from = Age, values_from = catch) %>% ungroup() %>% select(-year) %>% as.data.frame
  }
  if(type=='fyke'){
    CAA <- df %>% group_by(year) %>% summarize(catch_0=sum(catch_0),
                                               catch_1=sum(catch_1),
                                               catch_2=sum(catch_2),
                                               catch_3=sum(catch_3),
                                               catch_4=sum(catch_4),
                                               catch_5=sum(catch_5)) %>% select(-year)
  }
  return(all(CAA > 0))
}

fancy_scientific <- function(l) {
  if(max(l, na.rm=T) < 100){
    l <- format(l, scientific = FALSE, digits=2)
  } else {
    l <- format(l, scientific = TRUE, digits=2)
    l <- gsub("0e\\+00","0",l)
    l <- gsub("^(.*)e", "'\\1'e", l)
    l <- gsub("e\\+","e",l)
    l <- gsub("e", "%*%10^", l)
    l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  }
  parse(text=l)
}

plot_timeseries_subareaCV <- function(df, val='est', facet.var='Age', ytitle='Index of abundance (numbers)'){
  df[,'val'] <- df[,val]
  df[,'fvar'] <- factor(df[,facet.var], levels=unique(df[,facet.var]), labels=paste(facet.var,unique(df[,facet.var])))
  df$Year <- as.integer(df$year)

  df2 <- df %>% group_by(fvar, Subarea) %>% 
    summarize(meanCV = mean(se)) %>% 
    ungroup() %>% group_by(fvar, Subarea) %>% 
    summarize(labs = paste(format(round(meanCV,2), nsmall=2), collapse = "\n"), xpos=Inf, ypos=Inf) %>% as.data.frame
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=val)) + 
    ggplot2::ylab(ytitle) +
    ggplot2::geom_pointrange(ggplot2::aes(x=Year, ymin=lwr, ymax=upr, color=Subarea)) +
    ggplot2::geom_text(data=df2[df2$Subarea=='North',], aes(label=labs, x=xpos, y=ypos, color=Subarea), hjust = 1.2, vjust = 1.1, inherit.aes = FALSE) +
    ggplot2::geom_text(data=df2[df2$Subarea=='South',], aes(label=labs, x=xpos, y=ypos, color=Subarea), hjust = 1.2, vjust = 2.2, inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2])) +
    # ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(fvar), scales="free_y", ncol=3, strip.position = "top") +
    ggplot2::scale_x_continuous(expand=c(0.02,0.02), breaks = seq(min(df$Year),max(df$Year),2)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
  return(g)
}

plot_timeseries <- function(df, val='index', facet.var='age', relative=TRUE, 
                                    ytitle='Index (relative)', yzero=TRUE, addCI=FALSE, addCV=FALSE, alpha=0.05){
  df[,'val'] <- df[,val]
  if(relative) df[,'val'] <- df %>% group_by((!!as.symbol(facet.var))) %>% mutate(rel.val = val/mean(val)) %>% pull(rel.val)
  df[,'fvar'] <- factor(df[,facet.var], levels=unique(df[,facet.var]), labels=paste(facet.var,unique(df[,facet.var])))
  df$Year <- as.integer(df$year)
  
  # if(addCI){
  #   df <- df %>% group_by(fvar, Year) %>%
  #     summarize(lo = quantile(val, alpha/2, na.rm=TRUE), 
  #               hi = quantile(val, 1-alpha/2, na.rm=TRUE),
  #               cv = sd(val, na.rm=TRUE)/mean(val, na.rm=TRUE),
  #               val = mean(val, na.rm=TRUE)) %>% as.data.frame
  # }
  if(addCV){
    # df <- df %>% group_by(fvar) %>% mutate(xpos = 2020, ypos = 0.95*max(hi))
    # df2 <- df %>% group_by(fvar, cvar) %>% 
    #   summarize(meanCV = mean(cv), xpos=mean(xpos), ypos=mean(ypos)) %>% 
    #   ungroup() %>% group_by(fvar) %>% 
    #   summarize(labs = paste(round(meanCV,2), collapse = "\n"), xpos=mean(xpos), ypos=mean(ypos),
    #             ) %>% as.data.frame
    df2 <- df %>% group_by(fvar) %>% 
      summarize(meanCV = mean(se)) %>% 
      ungroup() %>% group_by(fvar) %>% 
      summarize(labs = paste(format(round(meanCV,2), nsmall=2), collapse = "\n"), xpos=Inf, ypos=Inf) %>% as.data.frame
  }
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=val)) + ggplot2::ylab(ytitle)  
  if(addCI){
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(x=Year, ymin=lwr, ymax=upr), color=NA, alpha=.15)
  }
  if(addCV){
    g <- g + ggplot2::geom_text(data=df2, aes(label=labs, x=xpos, y=ypos), hjust = 1.2, vjust = 1.1, inherit.aes = FALSE)
  }
  if(yzero){
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2]))
  } else {
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = function(x) range(pretty(x)))
  }
  g <- g + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(fvar), scales="free_y", ncol=3, strip.position = "top") +
    ggplot2::scale_x_continuous(expand=c(0.02,0.02), breaks = seq(min(df$Year),max(df$Year),2)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
  return(g)
}


plot_timeseries_compare <- function(df, val='index', facet.var='age', color.var='model', relative=TRUE, 
                                    ytitle='Index (relative)', yzero=TRUE, addCI=FALSE, addCV=FALSE, alpha=0.05){
  df[,'val'] <- df[,val]
  if(relative) df[,'val'] <- df %>% group_by((!!as.symbol(color.var)), (!!as.symbol(facet.var))) %>% mutate(rel.val = val/mean(val)) %>% pull(rel.val)
  df[,'fvar'] <- factor(df[,facet.var], levels=unique(df[,facet.var]), labels=paste(facet.var,unique(df[,facet.var])))
  df[,'cvar'] <- factor(df[,color.var], levels=unique(df[,color.var]), labels=unique(df[,color.var]))
  df$Year <- as.integer(df$year)

  # if(addCI){
  #   df <- df %>% group_by(fvar, cvar, Year) %>%
  #                summarize(lo = quantile(val, alpha/2, na.rm=TRUE), 
  #                          hi = quantile(val, 1-alpha/2, na.rm=TRUE),
  #                          cv = sd(val, na.rm=TRUE)/mean(val, na.rm=TRUE),
  #                          val = mean(val, na.rm=TRUE)) %>% as.data.frame
  # }
  if(addCV){
    # df <- df %>% group_by(fvar) %>% mutate(xpos = 2020, ypos = 0.95*max(hi))
    # df2 <- df %>% group_by(fvar, cvar) %>% 
    #   summarize(meanCV = mean(cv), xpos=mean(xpos), ypos=mean(ypos)) %>% 
    #   ungroup() %>% group_by(fvar) %>% 
    #   summarize(labs = paste(round(meanCV,2), collapse = "\n"), xpos=mean(xpos), ypos=mean(ypos),
    #             ) %>% as.data.frame
    df2 <- df %>% group_by(fvar, cvar) %>% 
      summarize(meanCV = mean(se)) %>% 
      ungroup() %>% group_by(fvar) %>% 
      summarize(labs = paste(format(round(meanCV,2), nsmall=2), collapse = "\n"), xpos=Inf, ypos=Inf) %>% as.data.frame
  }
  
  n.mods <- length(unique(df$cvar)) - 1
  cols <- c(RColorBrewer::brewer.pal(n = n.mods, name = "Set1"), "black")
  ltys <- c(rep(1, n.mods), 2)
  
  g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=val, color=cvar, linetype=cvar)) + ggplot2::ylab(ytitle)  
  # if(any(plot.opts$ci)){
  #   df <- df %>% dplyr::group_by(var) %>% dplyr::mutate(y_max = 1.2*max(val)) %>% as.data.frame   # trim large CIs to zoom to 120% of max MLE
  #   ind <- which(df$hi > df$y_max)
  #   df$hi[ind] = df$y_max[ind]
  # }

  if(addCI){
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(x=Year, ymin=lwr, ymax=upr, fill=cvar), color=NA, alpha=.15) + 
      ggplot2::scale_fill_manual(values=cols, name='')
  }
  if(addCV){
    g <- g + ggplot2::geom_text(data=df2, aes(label=labs, x=xpos, y=ypos), hjust = 1.2, vjust = 1.1, inherit.aes = FALSE)
  }
  if(yzero){
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2]))
  } else {
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = function(x) range(pretty(x)))
  }
  g <- g + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(fvar), scales="free_y", ncol=3, strip.position = "top") +
    ggplot2::scale_x_continuous(expand=c(0.02,0.02), breaks = seq(min(df$Year),max(df$Year),2)) + # breaks=scales::breaks_extended(5)
    scale_colour_manual(values=cols, name='') + scale_linetype_manual(values=ltys, name='') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
  return(g)
}

# --------------------------------------------------------------
# consistency (cohort tracking)
# function make cohorts
makecohorts <- function(mat){
  NAmat <- matrix(NA, NCOL(mat), NCOL(mat))
  matcoh1 <- rbind(NAmat, mat, NAmat)
  nr <- NROW(matcoh1)
  nc <- NCOL(mat)
  matcoh <- matrix(NA, nrow=nr, ncol=nc)
  for (i in 1:nc) matcoh[1:(nr-i),i] <- matcoh1[i:(nr-1),i]
  return(matcoh)
}

# function to plot all ages by all ages of data by cohorts
# assumes matrix has ages as columns and that cohorts are in rows
plotcoh <- function(matcoh, mytitle="", mylabels=NA){
  origpar <- par(no.readonly = TRUE)
  nc <- NCOL(matcoh)
  my.cor <- cor(matcoh, use="pairwise.complete.obs")
  my.cor.round <- round(my.cor,2)
  par(mfcol=c(nc,nc))
  par(oma=c(0,0,3,1),mar=c(1,1,0,0))
  for (i in 1:nc)
  {
    for (j in nc:1)
    {
      if (i == j)
      {
        plot(1:10,1:10,type='n',axes=FALSE)
        text(5,5, mylabels[i],cex=1.4)
      }
      if (i < j)
      {
        if (!is.na(my.cor[i,j]))
        {
          plot(matcoh[,i],matcoh[,j],axes=FALSE) # make sure have some data to plot
          xx <- matcoh[,i]
          yy <- matcoh[,j]
          my.fit <- lm(yy~xx)
          if (!is.na(my.fit$coefficients[2])) abline(my.fit,col="red")
          xrng <- data.frame(xx = seq(min(xx,na.rm=T),max(xx,na.rm=T),length.out=100))
          zz <- predict(my.fit, xrng, interval="confidence") 
          lines(xrng[,1],zz[,2],col="blue")
          lines(xrng[,1],zz[,3],col="blue")
        }
        if (is.na(my.cor[i,j]))
        {  # if not data, just make empty box
          plot(1:10,1:10,type='n',axes=FALSE)
        }
        box()
      }
      if (i > j)
      {
        plot(1:10,1:10,type='n',axes=FALSE)
        txt <- format(my.cor.round[i,j], nsmall=2)
        text(5,5,txt,cex = 4 * abs(my.cor.round[i,j]) / strwidth(txt) + 0.4*strwidth(txt))
        box()
      }
    }
  }
  title(mytitle, outer=TRUE)
}

plot_index_at_age_consistency <- function(df, col='rel.logind', plot=TRUE, fname='iaa-consistency'){
  models <- unique(df$model)
  n.mods = length(models)
  agelabs <- paste('Age', unique(df$Age))
  df$year = as.numeric(df$year)

  index.corr <- vector('list',n.mods)
  for(m in 1:n.mods){
    IAA <- df %>% filter(model == models[m]) %>%
              select(Age, year, (!!as.symbol(col))) %>% 
              complete(year = min(year):max(year)) %>%
              pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
              select(-c(year,`NA`)) %>%
              as.matrix
    coh <- makecohorts(IAA)
    index.corr[[m]] <- cor(coh, use="pairwise.complete.obs")

    if(plot){
      png(here('plots',paste0(fname,'-',col,'-',m,'.png')), units='in', width=7, height=7, res=150)
      plotcoh(coh, mytitle=paste0(models[m]," predicted"), mylabels=agelabs)
      dev.off()      
    }
  }
  return(index.corr)
}

plot_index_at_age_consistency_sdmtmb <- function(df, col='est', plot=TRUE, fname='iaa-consistency'){
  agelabs <- paste('Age', unique(df$Age))
  df$year = as.numeric(df$year)
  
    IAA <- df %>%
      select(Age, year, (!!as.symbol(col))) %>% 
      complete(year = min(year):max(year)) %>%
      pivot_wider(names_from = Age, values_from=(!!as.symbol(col))) %>%
      select(-c(year,`NA`)) %>%
      as.matrix
    coh <- makecohorts(IAA)
    index.corr <- cor(coh, use="pairwise.complete.obs")
    
    if(plot){
      # png(here('plots',paste0(fname,'.png')), units='in', width=7, height=7, res=150)
      png(fname, units='in', width=7, height=7, res=150)
      plotcoh(coh, mylabels=agelabs)
      dev.off()      
    }
  return(index.corr)
}


plot_correlation <- function(cohcor, models=NULL, alpha=0.05){
  df.cohcor <- cohcor %>% pivot_longer(-Age,
                                       names_to = 'Model',
                                       values_to = 'cor')
  df.cohcor <- df.cohcor %>% group_by(Age, Model) %>% 
    summarize(mean = round(mean(cor, na.rm=TRUE), 2),
              lo = round(quantile(cor, alpha/2, na.rm=TRUE), 2), 
              hi = round(quantile(cor, 1-alpha/2, na.rm=TRUE), 2)) %>% as.data.frame
  if(is.null(models)) models <- unique(df.cohcor$Model)
  df.cohcor$Model <- factor(df.cohcor$Model, levels=models, labels=models)
  n.mods <- length(models)
  cols <- c(RColorBrewer::brewer.pal(n = n.mods-1, name = "Set1"), "black")
  
  g <- ggplot(df.cohcor, aes(x=Age, y=mean, ymin=lo, ymax=hi, color=Model)) +
    geom_hline(yintercept=0, linetype=2) +
    geom_pointrange(size=1, position = position_dodge(width = 0.5)) +
    xlab('') +
    ylab('Correlation') +
    scale_colour_manual(values=cols, name='') +
    theme_bw(base_size = 14) +
    theme(legend.position = 'top')
  return(g)
}

plot_cv <- function(df, val='index.stand', models=NULL){
  df[,'val'] <- df[,val]

  df2 <- df %>% group_by(Age, model, year) %>%
    summarize(cv = sd(val, na.rm=TRUE)/mean(val, na.rm=TRUE),
              mean = mean(val, na.rm=TRUE)) %>% as.data.frame
  df2$cv[df2$cv > 2] = NA
  
  # df2 <- df %>% group_by(age, model, year) %>%
  #   summarize(cv = sd(val)/mean(val),
  #             mean = mean(val)) %>% ungroup() %>%
  #   group_by(age, model) %>% 
  #   summarize(meanCV = mean(cv)) %>% as.data.frame
    # pivot_wider(names_from = age, values_from=(!!as.symbol(col))) %>% as.data.frame
  
  if(is.null(models)) models <- unique(df2$model)
  df2$Model <- factor(df2$model, levels=models, labels=models)
  df2$Age <- as.factor(df2$Age)
  n.mods <- length(models)
  cols <- c(RColorBrewer::brewer.pal(n = n.mods-1, name = "Set1"), "black")
  
  g <- ggplot(df2, aes(x=Age, y=cv, fill=Model, group=interaction(Model,age))) +
    # geom_hline(yintercept=0, linetype=2) +
    # geom_point(size=3, position = position_dodge(width = 0.5)) +
    geom_boxplot(position = position_dodge(width=.75), alpha=0.4, coef = 0) +
    geom_point(shape=21, position = position_dodge(width=.75), size=2) +
    xlab('Age') +
    ylab('CV') +
    scale_y_continuous(expand=c(0.01,0.01), limits = function(x) c(0,range(pretty(x))[2])) +
    scale_fill_manual(values=cols, name='') +
    theme_bw(base_size = 14) +
    theme(legend.position = 'top')
  return(g)
}

# -----------------------------------------------------------------------------
# plot cohorts function from Johanna

# mat is a matrix of survey indices with ages in columns and years in rows
# change number of columns (ages) and the year range

## Number by age in time series ----
plot.cohort.year <- function(mat, age=1:16, survey.year=2003:2021, start.coh = 2003, stop.coh=2017,col1=1,col2=1, addPlot=F, ylab="Log10 (abundance index)"){
  ## Define first cohort: start.co
  ## Does the matrix include year and age information? Add if not
  mat <- mat[ ,colnames(mat) %in% age]
  mat[mat == 1] <- NA
  if(!is.null(rownames(mat))) survey.year <- as.numeric(rownames(mat))
  mat <- mat[row.names(mat) %in% survey.year,]
  age.mat <- matrix(age,nrow=nrow(mat),ncol=ncol(mat), byrow=T)
  sur.mat <- matrix(survey.year,nrow=nrow(mat),ncol=ncol(mat), byrow=F)
  coh.mat <- sur.mat - age.mat ### 
  if(addPlot == F) plot(sur.mat, mat, type="n", xlab="Year", ylab = ylab)
  for(i.y in start.coh:stop.coh){
    lines(sur.mat[coh.mat== i.y], mat[coh.mat == i.y], col=i.y, lwd=2)  
    text(sur.mat[coh.mat== i.y], mat[coh.mat == i.y], age.mat[coh.mat == i.y], col=i.y)
  }
}



