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
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH CITY, REDUCE AND SAVE
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS FOR 
# THE REDUCED OVERALL CUMULATIVE EXPOSURE-RESPONSE ONLY

coef <- coef1 <- coef2 <- coefint <- matrix(NA,nrow(cities),length(varper)+2,
  dimnames=list(cities$city))
vcov <- vcov1 <- vcov2 <- vcovint <- vector("list",nrow(cities))
names(vcov) <- names(vcov1) <- names(vcov2) <- names(vcovint) <- cities$city

################################################################################
# RUN THE MODEL FOR EACH CITY

# LOOP
for(i in seq(nrow(cities))) {

  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE CROSS-BASIS
  argvar <- c(vararg,list(knots=if(is.null(varper)) NULL else 
    quantile(data$tmean,varper/100,na.rm=T)))
  cb <- crossbasis(data$tmean,lag=lag,argvar=argvar,
    arglag=list(knots=logknots(lag,lagnk)),group=data$year)
  #summary(cb)
  
  # RUN THE MODELS WITHOUT INTERACTION
  model <- glm(formula,family=quasipoisson,data=data,na.action="na.exclude")
 
  # DEFINE INTERACTION TERMS, STANDARDIZED APPROX WITHIN 0-1
  # -INT1: CENTERED AT MID JUNE (OR DECEMBER FOR AUSTRALIA)
  # -INT2: CENTERED AT MID SEPTEMBER ((OR MARCH FOR AUSTRALIA))
  int1 <- ((data$dos-15)/120)*cb
  int2 <- ((data$dos-107)/120)*cb
  
  # RUN THE MODELS WITH INTERACTION
  model1 <- glm(formula1,family=quasipoisson,data=data)
  model2 <- glm(formula2,family=quasipoisson,data=data)
  
  # PREDICTION AND REDUCTION TO OVERALL CUMULATIVE EXPOSURE-RESPONSE
  # - AVERAGE, EARLY SUMMER, LATE SUMMER, INTERACTION TERMS
  # NB: CENTERING NOT NEEDED AT THIS STAGE, AS IT DOES NOT CHANGE COEF-VCOV
  red <- crossreduce(cb,model)
  coef[i,] <- coef(red)
  vcov[[i]] <- vcov(red)

  red1 <- crossreduce(cb,model1)
  coef1[i,] <- coef(red1)
  vcov1[[i]] <- vcov(red1)

  red2 <- crossreduce(cb,model2)
  coef2[i,] <- coef(red2)
  vcov2[[i]] <- vcov(red2)

  redint <- crossreduce(int1,model1)
  coefint[i,] <- coef(redint)
  vcovint[[i]] <- vcov(redint)
}

#

