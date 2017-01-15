################################################################################
# Updated version of the R code for the analysis in:
#
#   "Changes in susceptibility to heat within the summer:
#     a multi-country analysis"
#   Antonio Gasparrini and collaborators
#   American Journal of Epidemiology - 2016
#   http://www.ag-myresearch.com/2016_gasparrini_aje.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2016_gasparrini_AJE_Rcodedata
################################################################################

################################################################################
# MAIN RESULTS
################################################################################

################################################################################
# DESCRIPTIVE (TABLE 1 IN THE MANUSCRIPT)

# NUMBER OF LOCATIONS
length(dlist)

# TOTAL NUMBER OF DEATHS
sum(sapply(dlist,function(x) sum(x$death,na.rm=T)))

# STUDY PERIOD
(period <- rowMeans(sapply(dlist,function(x) range(x$year))))

# TEMPERATURE DISTRIBUTION IN EARLY AND LATE SUMMER (FIRST AND LAST 2 MONTHS)
rowMeans(sapply(dlist,function(x) summary(subset(x,dos<=61)$tmean)))
rowMeans(sapply(dlist,function(x) summary(subset(x,dos>61)$tmean)))

################################################################################
# EFFECTS BY COUNTRY (TABLE 2 IN THE MANUSCRIPT)
# NB: NOT IDENTICAL TO MANUSCRIPT, AS BASED ON UK ONLY

# FUNCTION FOR MULTIVARIATE WALD TEST
Ztest <- function(b1,V1,b2=NULL,V2=NULL) {
  invVp <- if(is.null(b2)) solve(V1) else solve(V1+V2)
  b <- if(is.null(b2)) b1 else b1-b2
  stat <- t(b)%*%invVp%*%(b)
  df <- length(b1)
  pchisq(stat,df,lower.tail=FALSE)
}

# MINIMUM MORTALITY PERCENTILE (MMP)
cenpercountry

# RR AT 90TH AND 99TH VS MMP (WITH 95%CI)
  # AVERAGE, EARLY AND LATE SUMMER
cp$allRRfit[predper==90];cp$allRRlow[predper==90];cp$allRRhigh[predper==90]
cp1$allRRfit[predper==90];cp1$allRRlow[predper==90];cp1$allRRhigh[predper==90]
cp2$allRRfit[predper==90];cp2$allRRlow[predper==90];cp2$allRRhigh[predper==90]

cp$allRRfit[predper==99];cp$allRRlow[predper==99];cp$allRRhigh[predper==99]
cp1$allRRfit[predper==99];cp1$allRRlow[predper==99];cp1$allRRhigh[predper==99]
cp2$allRRfit[predper==99];cp2$allRRlow[predper==99];cp2$allRRhigh[predper==99]

# MULTIVARIATE TESTS FOR A NULL INTERACTION (P-VALUE)
Ztest(coef(cpint),vcov(cpint))
Ztest(coef(cp1),vcov(cp1),coef(cp2),vcov(cp2))

#

