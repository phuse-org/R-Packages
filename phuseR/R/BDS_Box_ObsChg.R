
#' BDS_Box_ObsChg
#'
#' Combine outputs from BDS_Box_Obs and BDS_Box_Chg into a single page
#'
#' @export
#' @param bdsdset Dataset in the form of an R Data Frame
#' @param treatmentname Which treatment arm variable? e.g."TRTA" TRTA, TRTP, etc
#' @param useshortnames  #Rename Treatment Arms? (if you wish to display shorter names).TRUE OR FALSE
#' @param oldnames Treatment Arms old names .e.g "Xanomeline Low Dose","Xanomeline High Dose"
#' @param newnames Treatment Arms new names e.g., "X-low", "X-High"
#' @param usepopflag subset on a population flag. TRUE OR FALSE
#' @param popflag value "SAFFL"
#' @param testname test or parameter to be analyzed e.g."DIABP"
#' @param yaxislabel labels for y axis
#' @param selectedvisits visit numbers to be analyzed e.g 0,2,4,6,8,12,16,20,24
#' @param perpage how many visits to display per page
#' @param dignum number of digits in table, standard deviation = dignum +1
#' @param redoutliers True or False.
#' @param horizontallines True or False.
#' @param enterlimits True or False.
#' @param ANRLO lower limit(s)
#' @param ANRHI upper limit(s)
#' @param outputdirectory set output file directory
#' @param filetype output file type - TIFF or JPEG or PNG
#' @param pixelwidth choose output file size: pixel width
#' @param pixelheight choose output file size: pixel height
#' @param outputfontsize choose output font size
#' @param charttitle Title for the chart
#' @return PhUSE Figure 7.3 Box plot - Measurements and Change from Baseline by Analysis Timepoint, Visit and Treatment
#'
#' @import Hmisc
#' @import ggplot2
#' @import tools
#' @import gridExtra
#' @import data.table
BDS__Box_ObsChg<-function(bdsdset, treatmentname, useshortnames = c(TRUE,FALSE), oldnames, newnames,usepopflag = c(TRUE,FALSE), popflag, testname, yaxislabel, selectedvisits, perpage, dignum, redoutliers = c(TRUE, FALSE), horizontallines = c(TRUE,FALSE), enterlimits= c(TRUE,FALSE), ANRLO, ANRHI, outputdirectory, filetype = c("PNG","TIFF","JPEG"),  pixelwidth, pixelheight, outputfontsize, charttitle){


testresultsread <- bdsdset

#buildtable function to be called later, summarize data to enable creation of accompanying datatable
buildtable <- function(avalue, dfname, by1, by2, dignum){
  byvarslist <- c(by1,by2)
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
visitsplits <- 2*(ceiling((length(selectedvisits)/perpage)))
#for each needed page, subset selected visits by number per page
for(i in 1:visitsplits) {

  #subset on test, visits, population to be analyzed
  if (usepopflag == TRUE){
    testresults <- subset(testresultsread, PARAMCD == testname & AVISITN %in% selectedvisits[(initial):
                                                                                               (ifelse(perpage*i>length(selectedvisits),length(selectedvisits),perpage*i))]
                          & FLAG == "Y")

  } else {
    testresults <- subset(testresultsread, PARAMCD == testname & AVISITN %in% selectedvisits[(initial):(perpage*i)])
  }
  initial <- initial + perpage
  testresults<- data.table(testresults)

  #setkey for speed gains when summarizing
  setkey(testresults, USUBJID, TREATMENT, AVISITN)

  #CREATING THE 7.1 PLOT

  #create a variable for the out of limits data
  if (enterlimits == TRUE){
    testresults$OUT <- ifelse(testresults$AVAL < ANRLO | testresults$AVAL > ANRHI, testresults$AVAL, NA)
  } else if (enterlimits == FALSE){
    testresults$OUT <- ifelse(testresults$AVAL < testresults$ANRLO | testresults$AVAL > testresults$ANRHI, testresults$AVAL, NA)
  } else {print("WARNING - Manual entry of limits or automatic usage of limits in data not defined")}
  #specify plot
  p <- ggplot(testresults, aes(factor(AVISITN), fill = TREATMENT, AVAL))+ggtitle('7.1')
  # add notch, axis labels, legend, text size
  p1 <- p + geom_boxplot(notch = TRUE) + xlab("Visit Number") + ylab(yaxislabel) + theme(legend.position="bottom", legend.title=element_blank(),
                                                                                         text = element_text(size = outputfontsize),
                                                                                         axis.text.x  = element_text(size=outputfontsize),
                                                                                         axis.text.y = element_text(size=outputfontsize)) +ggtitle(charttitle)
  # add mean points
  p2 <- p1 + stat_summary(fun.y=mean, colour="dark red", geom="point", position=position_dodge(width=0.75))
  # out of limits jittered red points
  p3 <- p2 + geom_jitter(data = testresults, aes(factor(AVISITN), testresults$OUT), colour = "dark red", position = position_dodge(width=0.75))
  # horizontal limit lines
  if(enterlimits == TRUE){
    p4 <- p2 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
    pall <- p3 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
  } else if (enterlimits == FALSE) {
    p4 <- p2 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red")
    pall <- p3 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red")
  }
  #call summary table function
  summary <- buildtable(avalue = quote(AVAL), dfname= quote(testresults), by1 = "AVISITN", by2 = "TREATMENT", dignum)[order(AVISITN, TREATMENT)]
  table_summary <- data.frame(t(summary))

  t1theme <- ttheme_default(core = list(fg_params = list (fontsize = outputfontsize)))
  t1 <- tableGrob(table_summary, theme = t1theme, cols = NULL)

  # END 7.1 CREATION

  #CREATE 7.2

  chgp <- ggplot(testresults, aes(factor(AVISITN), fill = TREATMENT, CHG))+ggtitle('7.2')
  # add notch, axis labels, legend, text size
  chgp1 <- chgp + geom_boxplot(notch = TRUE) + xlab("Visit Number") + ylab(yaxislabel) + theme(legend.position="bottom", legend.title=element_blank(),
                                                                                         text = element_text(size = outputfontsize),
                                                                                         axis.text.x  = element_text(size=outputfontsize),
                                                                                         axis.text.y = element_text(size=outputfontsize)) +ggtitle(charttitle)
  # add mean points
  chgp2 <- chgp1 + stat_summary(fun.y=mean, colour="dark red", geom="point", position=position_dodge(width=0.75))

  # horizontal line at 0
  chgp3 <- chgp2 + geom_hline(yintercept = 0, colour = "red")

  if(enterlimits == TRUE){
    chgp4 <- chgp2 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
    chgpall <- chgp3 + geom_hline(yintercept = c(ANRLO,ANRHI), colour = "red")
  } else if (enterlimits == FALSE) {
    chgp4 <- chgp2 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red")
    chgpall <- chgp3 + geom_hline(yintercept = c(testresults$ANRLO,testresults$ANRHI), colour = "red")
  }

  #call summary table function
  chgsummary <- buildtable(avalue = quote(CHG), dfname= quote(testresults), by1 = "AVISITN", by2 = "TREATMENT", dignum)[order(AVISITN, TREATMENT)]
  chgtable_summary <- data.frame(t(chgsummary))

  chgt1theme <- ttheme_default(core = list(fg_params = list (fontsize = outputfontsize)))
  chgt1 <- tableGrob(chgtable_summary, theme = t1theme, cols = NULL)

  #create the plots

  #"U:","WINDOWS"

  if (filetype == "TIFF"){
    #Output to TIFF
    tiff(file.path(outputdirectory,paste("plot",i,".TIFF",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) {
      grid.arrange(pall, chgpall, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, chgp3, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, chgp4, t1,chgt1, ncol = 2, nrow = 2)
    } else {
      grid.arrange(p2, chgp2, t1,chgt1, ncol = 2, nrow = 2)
    }
    dev.off()
  }
  if (filetype == "JPEG") {
    # Optionally, use JPEG
    jpeg(file.path(outputdirectory,paste("plot",i,".JPEG",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) {
      grid.arrange(pall,chgpall, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, chgp3, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, chgp4, t1,chgt1, ncol = 2, nrow = 2)
    } else {
      grid.arrange(p2, chgp2, t1,chgt1, ncol = 2, nrow = 2)
    }
    dev.off()
  }
  if (filetype == "PNG") {
    # Optionally, use PNG
    png(file.path(outputdirectory,paste("plot",i,".PNG",sep = "" )), width = pixelwidth, height = pixelheight, units = "px", pointsize = 12)
    if (redoutliers == TRUE & horizontallines == TRUE) {
      grid.arrange(pall, chgpall, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == TRUE & horizontallines == FALSE) {
      grid.arrange(p3, chgp3, t1,chgt1, ncol = 2, nrow = 2)
    } else if (redoutliers == FALSE & horizontallines == TRUE) {
      grid.arrange(p4, chgp4, t1,chgt1, ncol = 2, nrow = 2)
    } else {
      grid.arrange(p2, tchgp2, t1,chgt1, ncol = 2, nrow = 2)
    }
    dev.off()
  }
}
}

