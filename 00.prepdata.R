################################################################################
# Updated version of the R code for the analysis in:
#
#   "Changes in susceptibility to heat within the summer:
#     a multi-country analysis"
#   Antonio Gasparrini and collaborators
#   American Journal of Epidemiology - 2016
#
# Update: 09 May 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
#
# * an updated version of this code, (hopefully) compatible with future
#   versions of the software, is available at the personal website of the
#   first author (www.ag-myresearch.com)
################################################################################

# This code reproduces the analysis with the subset of data only including UK
# For this reason, the results are slightly different from those published

################################################################################
# PREPARE THE DATA
################################################################################

# LOAD THE PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines)

# CHECK VERSION OF THE PACKAGE
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")

# LOAD THE DATASET (INCLUDING THE 10 UK REGIONS ONLY)
orig <- read.csv("regEngWales19902012.csv")
orig$date <- as.Date(orig$date)

# ARRANGE THE DATA AS A LIST OF DATA SETS
regions <- as.character(unique(orig$regnames))
dlist <- lapply(regions,function(x) orig[orig$regnames==x,])
names(dlist) <- regions

# METADATA FOR LOCATIONS
cities <- data.frame(
  city = regions,
  cityname = c("North East","North West","Yorkshire & Humber","East Midlands",
    "West Midlands","East","London","South East","South West","Wales"),
  lat = c(55.00000,54.01670,53.56670,52.98000,52.47000,52.24000,
    51.50720,51.31670,50.96000,52.32044),
  long = c(-1.866700,-2.633300,-1.200000,-0.750000,-2.290000,0.410000,-0.127500,
    -0.500000,-3.220000,-3.732616)
)

# IDENTIFY THE 4 AVERAGE HOTTEST MONTHS IN EACH LOCATION
t(sapply(dlist,function(x) 
  sort((1:12)[order(tapply(x$tmean,x$month,mean,na.rm=T))][9:12])))

# RESTRICT TO SUMMER AND CREATE DAY OF SEASON VARIABLE
dlist <- lapply(seq(dlist),function(i) {
  sub <- subset(dlist[[i]],month %in% 6:9)
  sub$dos <- sequence(tapply(sub$date,sub$year,length))
  sub
})

################################################################################

# PARAMETERS FOR THE EXPOSURE-RESPONSE FUNCTION
vararg <- list(fun="bs",degree=2)
varper <- c(75)

# PARAMETERS FOR THE LAG-RESPONSE FUNCTION
lag <- 10
lagnk <- 2
lagdf <- 4

# DEGREES OF FREEDOM FOR SEASONALITY
dfseas <- 4

# DEGREES OF FREEDOM/DECADE FOR TIME TREND
dftrend <- 1

# MODEL FORMULAE
formula <- death ~ cb + dow + ns(dos,df=dfseas):factor(year) +
  ns(date,df=round(length(unique(year))/dftrend/10))
formula1 <- death ~ cb + dow + ns(yday,df=dfseas):factor(year) +
  ns(date,df=round(length(unique(year))/dftrend/10)) + int1
formula2 <- death ~ cb + dow + ns(yday,df=dfseas):factor(year) +
  ns(date,df=round(length(unique(year))/dftrend/10)) + int2

#
