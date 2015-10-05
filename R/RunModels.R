#' RunModels   
#' 
#' Run cross-validation or bootstrap in series or parallel 
#' @name  RunModels
#' 
#' @param form.list named list of formulas  
#' @param data data matrix  
#' @param classifier character list of classification models. See output of 
#' \code{names(TrainModels())}. 
#' @param cv,nBoots number of cv folds, bootsraps 
#' @param para named list of model hyper-parameters
#' @param opt.para (logical) perform parameter optimization ?
#' @param val.propo proportion of test set to set aside for seelection of classification threshold
#' @param return.model (logical) return trained model and predicted probabilities? 
#' @param varimp (logical) return variable importance ?
#' @param parallel (logical) train and test models in parallel using parallel package ? 
#' @param mc.cores number of cores to use in parallel mode 
#' @param seed random seed for reproducibility 
#' @param \dots further arguments passed to or from other methods.
#' @return list of length cv or nBoots containing performance measures for each 
#' classifier as a named list with  
#' \item{model}{trained model}
#' \item{para}{named list of model hyper-paramters (tunned values if opt.para = TRUE)}
#' \item{run.time}{compute time} 
#' \item{varimp}{variable importance if classifier is GLM, GBM or RF}
#' \item{perf}{a list with two data.frames: val.perf and tst.perf containing performance measures 
#' for validation and test sets }
#' 
#'  @details
#' \enumerate{
#'   \item \code{RunCV}:  
#'   This function runs in cross-validation style to train, test and select optimal 
#'   model hyper-parameters. A grid search is performed for each cross-validation using the grid 
#'   paramerter in \code{para} and the best parameters used to train the final model.
#'   \item \code{RunBootstrap}: Runs models in bootstrap (balanced) style to train, test and if desired also 
#'   perform paramter optimization for each bootstrap just as in \code{RunCV}. \code{RunCV} can be run first and 
#'   pass the output to \code{get.para} to select the para returned by the cross-validation with the 
#'   best performance. 
#' }
#'

#' @author Che Ngufor <Ngufor.Che@@mayo.edu>
#'
#' @import ELR caret
#'
NULL 
#' @rdname RunModels   
#' @export
#' @examples
#' \dontrun{
#' para = list(
#' method = "cv", 
#' tuneLength = 4,
#' number = 3,
#' GLM = c(), 
#' ELR = list(grid = expand.grid(gamma = c(0.01, 0.25), p =  c(100, 400)),  
#'            ken = "sigmoid", p = 400, gamma = 10.05), 
#' GBM = list(n.trees=10, interaction.depth=1, shrinkage=0.01, 
#'            mstop = 100, n.minobsinnode = 10,prune = "no", 
#'            grid = expand.grid(interaction.depth = c(1, 3, 5), 
#'                    n.trees= c(10, 20, 50), shrinkage=c(0.01), 
#'                    n.minobsinnode=10)),
#' RF = list(ntree = 100, mtry = 5), 
#' pLDA = list(lambda = 0.5),
#' avNNET = list(nnet.size = 50, nnet.decay= 0.001, nnet.bag = 10, maxit = 200), 
#' SVM = list(C = .05, sigma = 0.5, kernel= "rbfdot", grid = expand.grid(C = c(0.01, 1, 10, 100), 
#'                                  sigma = c(0.001, 0.05, 0.5))),
#' GLMnet = list(alpha = 1, lambda = 1, grid = expand.grid(alpha=c(0,0.5, 1),
#'              lambda= c(0.01, 0.05, 0.1, 0.5)))
#' ) 
#' # Pima Indian diabetes data 
#' link <-"http://archive.ics.uci.edu/ml/machine-learning-databases/
#' pima-indians-diabetes/pima-indians-diabetes.data"
#' dat <- read.table(file = link, header = FALSE, sep = ",") 
#' rhs.vars <- c("n.preg","glu","bp","skin","insulin","bmi","ped","age") ## independent variables 
#' outcome <- c("class1", "class2") ## outcome: diabetic or not
#' dat <- cbind(dat, dat[, 9]) 
#' names(dat) <- c(rhs.vars, outcome) 
#' dat[, rhs.vars] <-  normalize01(dat[, rhs.vars])
#' 
#' form <- lapply(outcome, function(y) as.formula(paste0(paste0(y, "~"), 
#' paste0(rhs.vars, collapse= "+"))))
#' names(form) <- outcome
#' classifier <- c("GLM", "GBM" , "ELR", "avNNET", "GLMnet", "SVM", "RF")
#' 
#' cv.res <- RunCV(form=form, data=dat, classifier=classifier, cv = 5, 
#'                 para=para, opt.para=TRUE, val.propo = 0.10, return.model = TRUE, 
 #'                varimp = TRUE, parallel = FALSE, mc.cores = cores, seed=12345)
#' cv.para <- get.para(object=cv.res, para=para)
#' cv.perf <-  get.cvResults(object=cv.res)
#' 
#' boot.res <- RunBootstrap(form=form, data=dat, classifier=classifier,  
#'                   nBoots = 5, para=cv.para, opt.para=FALSE, val.propo = 0.10, 
#'                   return.model = FALSE, varimp = FALSE, parallel = FALSE, 
#'                     mc.cores = cores, seed=12345)
#' boot.perf <-  get.bootResults(object=boot.res)
#' 
#'}
### cross-validation 
RunCV <- function(form.list, data, classifier, cv = 2, para, opt.para=FALSE, 
          val.propo = 0.10, return.model = FALSE, varimp = TRUE, parallel = TRUE, 
          mc.cores = getOption("mc.cores", 2L), seed=12345, ...){

resp.vars <- names(form.list)

XX <- data.frame(model.matrix(object = form.list[[1]], data = data)[, -1, drop = FALSE])
rhs.vars <- colnames(XX) 
cons.vars <- sapply(XX, Cons.Vars, thresh = 0.001) ### drop vary rare binary variables
drop.vars <- rhs.vars[cons.vars]
rhs.vars <- rhs.vars[!cons.vars]
XX <- XX[, rhs.vars]
data <- cbind.data.frame(data[, resp.vars, drop=FALSE], XX)


if(parallel) {
set.seed(seed, "L'Ecuyer")          
ix.cv <- Split(Nt = dim(data)[1], K= cv)
pfun <-  get("mclapply")
} else {
set.seed(seed) 
ix.cv <- Split(Nt = dim(data)[1], K= cv)
pfun = get("lapply")
}

CV.res <- pfun(1:cv, function(kk, ...){ 
 
dd.dat <- SplitCrossVal(data, kk, ix.cv, val.propo = val.propo)
dat.trn  <- dd.dat$dat.trn 
dat.val <-  dd.dat$dat.val 
dat.tst <- dd.dat$dat.tst 

resp.res <-  tryCatch(
{

lapply(resp.vars, function(yy){
### classification 
Y.trn <- dat.trn[, yy, drop = FALSE]
X.trn <-  dat.trn[, rhs.vars,drop=FALSE]
Y.val <- dat.val[, yy, drop = FALSE]
X.val <-  dat.val[, c(yy, rhs.vars),drop=FALSE]
Y.tst <- dat.tst[, yy, drop = FALSE]
X.tst <-  dat.tst[, rhs.vars,drop=FALSE]

cls.res <- lapply(classifier, function(x) Train.Validate.Test(classifier=x, X.trn=X.trn, Y.trn=Y.trn,
              X.val=X.val,Y.val=Y.val,X.tst=X.tst,Y.tst=Y.tst,para=para,opt.para=opt.para, 
              return.model=return.model,varimp=varimp))      
collect.garbage() 
names(cls.res) <- classifier                 
cls.res
}
)
} , error=function(e){ 
cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
": original error message = ", e$message, "\n") 
list()
}
) ## tryCatch

cat("Done CV :", kk,  "\n")
   
names(resp.res)<- resp.vars 
resp.res
}, mc.cores = mc.cores)  ## pfun  

return(list(CV.res = CV.res, drop.vars = drop.vars, rhs.vars = rhs.vars, 
resp.vars = resp.vars, return.model =return.model, varimp=varimp))
}

#'
#' @rdname RunModels   
#' @export
# bootstrap
RunBootstrap <- function(form.list, data, classifier, nBoots = 5, para, opt.para=FALSE, 
             val.propo = 0.10, return.model = FALSE, varimp = FALSE, parallel = TRUE, 
             mc.cores = getOption("mc.cores", 2L),  seed=12345, ...){

  resp.vars <- names(form.list)
  XX <- data.frame(model.matrix(object = form.list[[1]], data = data)[, -1, drop = FALSE])
  rhs.vars <- colnames(XX) 
  cons.vars <- sapply(XX, Cons.Vars, thresh = 0.001)  ### drop vary rare binary variables 
  drop.vars <- rhs.vars[cons.vars]
  rhs.vars <- rhs.vars[!cons.vars]
  XX <- XX[, rhs.vars]
  data <- cbind.data.frame(data[, resp.vars, drop=FALSE], XX)
    
nobs<-nrow(data)
rownames(data) <- NULL 
ix.boot <- 1:nobs

if(parallel) {
set.seed(seed, "L'Ecuyer")   
idx <- sample(matrix(t(replicate(nBoots, ix.boot))))      
pfun <-  get("mclapply")
} else {
set.seed(seed) 
idx <- sample(matrix(t(replicate(nBoots, ix.boot))))
pfun = get("lapply")
}


Boot.res <- pfun(1:nBoots, function(kk, ...){ 

inbag = idx[((kk-1)*nobs+1):(kk*nobs)]
outbag <- setdiff(ix.boot, inbag)

dat.trn <- data[inbag, ,drop = FALSE]
XX <- data[outbag, ,drop = FALSE] 
ix <- sample(nrow(XX), floor(nrow(XX)*val.propo) )
dat.val <- XX[ix, ,drop = FALSE] 
dat.tst <- XX[-ix, ,drop = FALSE] 
 
resp.res <-  tryCatch(
{

lapply(resp.vars, function(yy){

if(yy%in%names(para) ) para <- para[[yy]]

### classification 
Y.trn <- dat.trn[, yy, drop = FALSE]
X.trn <-  dat.trn[, rhs.vars,drop=FALSE]
Y.val <- dat.val[, yy, drop = FALSE]
X.val <-  dat.val[, c(yy, rhs.vars),drop=FALSE]
Y.tst <- dat.tst[, yy, drop = FALSE]
X.tst <-  dat.tst[, rhs.vars,drop=FALSE]

cls.res <- lapply(classifier, function(x) Train.Validate.Test(classifier=x, X.trn=X.trn, Y.trn=Y.trn,
              X.val=X.val,Y.val=Y.val,X.tst=X.tst,Y.tst=Y.tst,para=para,opt.para=opt.para, 
              return.model=return.model,varimp=varimp))      
collect.garbage() 
names(cls.res) <- classifier                 
cls.res
}
)
} , error=function(e){ 
cat("Error in the Expression: ",  paste(e$call, collapse= ", "), 
": original error message = ", e$message, "\n") 
list()
}
) ## tryCatch

cat("Done bootstrap :", kk,  "\n")  
names(resp.res)<- resp.vars 
resp.res
}, mc.cores = mc.cores)  ## pfun  

return(list(Boot.res = Boot.res, drop.vars = drop.vars, rhs.vars = rhs.vars, 
resp.vars = resp.vars, return.model =return.model, varimp=varimp))
}








