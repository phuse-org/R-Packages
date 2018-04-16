
#' BDS_Box_LMM
#'
#' Last/Min/Max Baseline and Last/Min/Max Postbaseline Measures across multiple studies
#'
#' @export
#' @param bdsdset Dataset in the form of an R Data Frame
#' @param treatmentname Which treatment arm variable? e.g."TRTA" TRTA, TRTP, etc
#' @param useshortnames  #Rename Treatment Arms? (if you wish to display shorter names).TRUE OR FALSE
#' @param oldnames Treatment Arms old names .e.g "Xanomeline Low Dose","Xanomeline High Dose"
#' @param newnames Treatment Arms new names e.g., "X-low", "X-High"
#' @param usepopflag subset on a population flag. TRUE OR FALSE
#' @param popflag value "SAFFL"
#' @param crit MIN, MAX, or LAST
#' @param testname test or parameter to be analyzed e.g."DIABP"
#' @param yaxislabel labels for y axis
#' @param selectedstudies Study names to be analyzed e.g 0,2,4,6,8,12,16,20,24
#' @param perpage how many visits to display per page
#' @param dignum number of digits in table, standard deviation = dignum +1
#' @param redoutliers True or False.
#' @param horizontallines True or False.
#' @param enterlimits True or False.
#' @param ANRLO lower limit(s)
#' @param ANRHI upper limit(s)
#' @param outputdirectory set output file directory
#' @param filetype output file type - TIFF or JPEG or PNG
#' @param pixelwidth choose output file size: pixel width (eg 1200)
#' @param pixelheight choose output file size: pixel height (eg 1000)
#' @param outputfontsize choose output font size
#' @param charttitle Title for the chart
#' @return PhUSE Figure 7.6 Box Plot - Last/Min/Max Baseline and Last/Min/Max Postbaseline Measures
#'
#' @import Hmisc
#' @import ggplot2
#' @import tools
#' @import gridExtra
#' @import data.table
#' @import dplyr
BDS__Box_LMM<-function(bdsdset, treatmentname, useshortnames = c(TRUE,FALSE), oldnames, newnames,usepopflag = c(TRUE,FALSE), popflag, crit, testname, yaxislabel, selectedstudies, perpage, dignum, redoutliers = c(TRUE, FALSE), horizontallines = c(TRUE,FALSE), enterlimits= c(TRUE,FALSE), ANRLO, ANRHI, outputdirectory, filetype = c("PNG","TIFF","JPEG"),  pixelwidth, pixelheight, outputfontsize, charttitle){
  
testresultsread <- bdsdset

#buildtable function to be called later, summarize data to enable creation of accompanying datatable
buildtable <- function(avalue, dfname, by1, by2, by3, dignum){
  byvarslist <- c(by1,by2,by3)
  summary <- eval(dfname)[,list(
    n = .N,
    mean = round(mean(eval(avalue), na.rm = TRUE), digits=dignum),
    sd = round(sd(eval(avalue), na.rm = TRUE), digits=dignum+1),
    min = round(min(eval(avalue), na.rm = TRUE), digits=dignum),
    q1 = round(quantile(eval(avalue), .25, na.rm = TRUE), digits=dignum),
    mediam = round(median(eval(avalue), na.rm = TRUE), digits=dignum),
    q3 = round(quantile(eval(avalue), .75, na.rm = TRUE), digits = dignum),
    max = round(max(eval(avalue), na.rm = TRUE), digits = dignum)
  ), 
  by = byvarslist]
  
  return(summary)
}

#SELECT VARIABLES (examples in parenthesis): TREATMENT (TRTP, TRTA), PARAMCD (LBTESTCD)
#colnames(testresults)[names(testresults) == "OLD VARIABLE"] <- "NEW VARIABLE"

colnames(testresultsread)[names(testresultsread) == treatmentname] <- "TREATMENT"

colnames(testresultsread)[names(testresultsread) == popflag] <- "FLAG" #select population flag to subset on such as SAFFL or ITTFL

if (useshortnames == TRUE){
  for(i in 1:length(oldnames)) {
    testresultsread$TREATMENT <- ifelse(testresultsread$TREATMENT == oldnames[i], as.character(newnames[i]), as.character(testresultsread$TREATMENT))
  }
}

#determine number of pages needed
initial <- 1
studysplits <- ceiling((length(selectedstudies)/perpage))
#for each needed page, subset selected visits by number per page
for(i in 1:studysplits) {
  
  if ('Total' %in% selectedstudies){
    total <- testresultsread
    totald <-data.table(total)
    totald$STUDYID <- 'Total'
    testresultsread2 <- bind_rows(testresultsread,totald)
  }
  else{
    testresultsread2 <- testresultsread 
  }
  
  #subset on test, visits, population to be analyzed
  if (usepopflag == TRUE){
    testresults <- subset(testresultsread2, PARAMCD == testname & STUDYID %in% selectedstudies[(initial):
                                                                                               (ifelse(perpage*i>length(selectedstudies),length(selectedstudies),perpage*i))] 
                          & FLAG == "Y")
    
  } else {
    testresults <- subset(testresultsread2, PARAMCD == testname & STUDYID %in% selectedstudies[(initial):(perpage*i)]) 
    
  }
  if (usepopflag == TRUE){
    testresults <- subset(testresultsread2, PARAMCD == testname & STUDYID %in% selectedstudies[(initial):
                                                                                               (ifelse(perpage*i>length(selectedstudies),length(selectedstudies),perpage*i))] 
                          & FLAG == "Y")
    
  } else {
    testresults <- subset(testresultsread2, PARAMCD == testname & STUDYID %in% selectedstudies[(initial):(perpage*i)])  
  }
  initial <- initial + perpage
  testresults2<- data.table(testresults)
  
  if (crit == 'MAX'){
    
    pretrt1 <- filter(testresults2, testresults2$ABLFL == 'Y' & testresults2$ANL01FL == 'Y')
    
    pretrt <- group_by(pretrt1, USUBJID,STUDYID) %>%
      filter(AVAL == max(AVAL)) %>%
      filter( 1:n() == 1)
    
    
    
    posttrt1 <- filter(testresults2,testresults2$ADT >= TRTSDT & testresults2$ANL01FL == 'Y')
    
    posttrt <- group_by(posttrt1, USUBJID,STUDYID) %>%
      filter(AVAL == max(AVAL)) %>%
      filter( 1:n() == 1)
    
    
  }
  
  if (crit == 'MIN'){
    
    pretrt1 <- filter(testresults2, testresults2$ABLFL == 'Y' & testresults2$ANL01FL == 'Y')
    
    pretrt <- group_by(pretrt1, USUBJID,STUDYID) %>%
      filter(AVAL == min(AVAL)) %>%
      filter( 1:n() == 1)
    
    posttrt1 <- filter(testresults2,testresults2$ADT >= TRTSDT & testresults2$ANL01FL == 'Y')
    
    posttrt <- group_by(posttrt1, USUBJID,STUDYID) %>%
      filter(AVAL == min(AVAL)) %>%
      filter( 1:n() == 1)
  }
  
  if (crit == 'LAST'){
    
    pretrt1 <- filter(testresults2, testresults2$ABLFL == 'Y' & testresults2$ANL01FL == 'Y')
    
    pretrt <- group_by(pretrt1, USUBJID,STUDYID) %>%
      filter(ADT == max(ADT)) %>%
      filter( 1:n() == 1)
    
    posttrt1 <- filter(testresults2,testresults2$ADT >= TRTSDT & testresults2$ANL01FL == 'Y')
    
    posttrt <- group_by(posttrt1, USUBJID,STUDYID) %>%
      filter(ADT == max(ADT)) %>%
      filter( 1:n() == 1)
  }
  
  pretrt$PERIOD <- 'BASE'
  
  posttrt$PERIOD <- 'POST'
  
  testresults3 <- bind_rows(pretrt,posttrt)
  
  testresults <- data.table(testresults3)
  
  #setkey for speed gains when summarizing
  setkey(testresults, USUBJID, TREATMENT, STUDYID, PERIOD)
  
  #create a variable for the out of limits data
  if (enterlimits == TRUE){
    testresults$OUT <- ifelse(testresults$AVAL < ANRLO | testresults$AVAL > ANRHI, testresults$AVAL, NA)
  } else if (enterlimits == FALSE){
    testresults$OUT <- ifelse(testresults$AVAL < testresults$ANRLO | testresults$AVAL > testresults$ANRHI, testresults$AVAL, NA)
  } else {print("WARNING - Manual entry of limits or automatic usage of limits in data not defined")}
  #specify plot
  p <- ggplot(testresults, aes(factor(PERIOD), fill = TREATMENT, AVAL))
  # add notch, axis labels, legend, text size
  p1 <- p + geom_boxplot(notch = TRUE) + xlab("Period") + ylab(yaxislabel) + theme(legend.position="bottom", legend.title=element_blank(), 
                                                                                   text = element_text(size = outputfontsize),
                                                                                   axis.text.x  = element_text(size=outputfontsize),
                                                                                   axis.text.y = element_text(size=outputfontsize)) +ggtitle(charttitle)
  # add mean points
  p2 <- p1 + stat_summary(fun.y=mean, colour="dark red", geom="point", position=position_dodge(width=0.75))
  # out of limits jittered red points
  p3 <- p2 + geom_jitter(data = testresults, aes(factor(PERIOD), testresults$OUT), colour = "dark red", position = position_dodge(width=0.75))
  # horizontal limit lines
  if(enterlimits == TRUE){
    p4 <- p2 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
    pall <- p3 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
  } else if (enterlimits == FALSE) {
    p4 <- p2 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red")
    pall <- p3 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red") 
  }
  #call summary table function
  summary <- buildtable(avalue = quote(AVAL), dfname= quote(testresults), by1 = "STUDYID", by2 = "PERIOD", by3 = "TREATMENT", dignum)[order(STUDYID, PERIOD, TREATMENT)]
  table_summary <- data.frame(t(summary))           
  
  t1theme <- ttheme_default(core = list(fg_params = list (fontsize = outputfontsize)))
  t1 <- tableGrob(table_summary, theme = t1theme, cols = NULL) 
  
  if (filetype == "TIFF"){
    #Output to TIFF
    tiff(file.path(outputdirectory,paste("plot_6",i,".TIFF",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) { 
      grid.arrange(pall, t1, ncol = 1)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, t1, ncol = 1)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, t1, ncol = 1)
    } else {
      grid.arrange(p2, t1, ncol = 1)
    }
    dev.off()
  }
  if (filetype == "JPEG") { 
    # Optionally, use JPEG
    jpeg(file.path(outputdirectory,paste("plot_6",i,".JPEG",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) { 
      grid.arrange(pall, t1, ncol = 1)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, t1, ncol = 1)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, t1, ncol = 1)
    } else {
      grid.arrange(p2, t1, ncol = 1)
    }
    dev.off()
  }
  if (filetype == "PNG") { 
    # Optionally, use PNG
    png(file.path(outputdirectory,paste("plot_6",i,".PNG",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) { 
      grid.arrange(pall, t1, ncol = 1)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, t1, ncol = 1)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, t1, ncol = 1)
    } else {
      grid.arrange(p2, t1, ncol = 1)
    }
    dev.off()
  }
}
}