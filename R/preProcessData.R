#' preProcessData
#' 
#' Preprocess data: transform, normalize, select, impute variables etc.  
#'
#' @param file.name  data file name
#' @param impute (logical) impute missing values with missforest or drop missing values. 
#' Note imputation with \code{\link[missForest]{missForest}} is going to take a considerable amount of time.  
#' @return  training data, response variables, predictors and formula.  
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' 
#' @import missForest
#'
#' @export 
preProcessData <- function(file.name= "", impute = FALSE){
dat <- read.table(file = file.name, header = TRUE, sep = ",") 
dat <- dat[, -c(1,2)]
ix.vars <- c(1,3,4,8,9,10,14:24,26:87,152,155:162,167, 169, 170,171,172,173,174,220,221,
              222, 242,256, 267, 268, 269, 275, 281, 287, 293)
dt <- dat[-1,ix.vars]
### Discharge destination (home or non-Home) 
dt$Rehab <- ifelse(dt$dischdest%in%c("Facility Which was Home", "Home"), 0, 1)
dt$Readmission <- iv <- ifelse(dt$unplannedreadmission1 == "Yes" | dt$unplannedreadmission2 == "Yes" | 
                         dt$unplannedreadmission3 == "Yes" | dt$unplannedreadmission4 == "Yes" |
                         dt$unplannedreadmission5 == "Yes", 1, 0) 
resp.vars <- c("Rehab", "Readmission")
rhs.vars <- unique(c("race_new", "ethnicity_hispanic", "inout",  "transt", "age", "anesthes","attend", 
              "surgspec", "height", "weight","diabetes" ,"smoke","packs", "etoh","dyspnea", "dnr", 
              "fnstatus2", "ventilat", "hxcopd", "cpneumon","ascites","esovar","hxchf","hxmi", 
              "prvpci", "prvpcs","hxangina","hypermed","hxpvd","tumorcns", "para","quad","discancr",
              "wndinf","steroid", "wtloss", "bleeddis","dprbun", "dprcreat","dpralbum","dprbili" ,
              "dprsgot","dpralkph","dprbun", "dprcreat","dpralbum","dprbili" ,"dprsgot","dpralkph",
              "dprwbc","dprhct" ,"dprplate","dprptt","dprpt", "dprinr", "prsodm","prbun","prcreat",
              "pralbum","prbili","prsgot", "pralkph", "prwbc","prhct","prplate","prptt","prinr",
              "prpt", "asaclas"))

dt$age[dt$age == "90+"] <- 90
dt$age[dt$age == "-99"] <- NA
dt$age <- as.numeric(dt$age)

dt$age <- cut(dt$age, breaks = c(18, 25, 35, 45, 55, 65, 75, 85, 95), 
labels = c("18.to.24","25.to.34","35.to.44", "45.to.54", "55.to.64", "65.to.74", 
           "75.to.84", "85.plus"), include.lowest = FALSE,  right = FALSE,  
           ordered_result = FALSE)

clss <- sapply(dt[, rhs.vars], class)
ix.cls <- names(clss[(clss%in%c("integer", "numeric"))])
nas <- dt[, ix.cls] == -99
dt[, ix.cls][nas] <- NA 
### drop columns with more than 30% missing 
cont.vars <- ix.cls[!sapply(dt[, ix.cls], function(x) (sum(is.na(x))/nrow(dt)) > 0.3)] 
### normalize continous variables to (0, 1) 
dt[, cont.vars] <- normalize01(dt[, cont.vars])

#### characters 
ix.cls2 <- names(clss[(clss%in%c("factor","character"))])
dt[, ix.cls2] <- do.call(cbind.data.frame, lapply(dt[, ix.cls2], as.factor))

dt <- dt[, c(resp.vars, cont.vars, ix.cls2)]  
rhs.vars <- c(cont.vars,ix.cls2) 
form.list <- lapply(resp.vars, function(x) as.formula(paste(paste(x, "~"), paste(rhs.vars, collapse= "+"))))                        
names(form.list) <- resp.vars

if(impute){
Ximp <- missForest(dt, verbose = TRUE)
XY.dat <- Ximp$ximp
} else XY.dat <- na.omit(dt)

list(XY.dat = XY.dat, resp.vars = resp.vars, rhs.vars = rhs.vars, form.list = form.list)
}



