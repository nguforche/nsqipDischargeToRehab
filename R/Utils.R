#' Utility functions.
#' 
#' Some useful miscellaneous functions 
#'
#' @name Utils
#' @param D.dat,x A matix/data frame 
#' @param Nt Sample size 
#' @param K,k number of cross-validation and cross-validation number 
#' @param ix.cv cross-validation indices 
#' @param val.propo training set proportion 
#' @param dummy  numeric or factor variable 
#' @param thresh threshold for deciding if variable is constant 
#' @param threshold,prevalence classification threshold and observed class prevalence
#' @param prob,pred,obs predicted probabilities and true observed binary response vectors 
#' @param uniqueCut minimum number of unique values for a variable to not be 
#' considered as constant
#' @param object train cv or bootstrap models
#' @param para model hyper-parameters 
#' @param alpha significance level 
#' @param top top number of variables to return 
#' @param deparse.level see \code{\link[base]{cbind}} 
#' @param row.names,check.rows,check.names,stringsAsFactors = FALSE see \code{\link[base]{data.frame}} 
#' @param \dots further arguments passed to or from other methods.
#' @return  training, validation, test data sets or K-fold cross-validation data sets. 
#'
#'  @details
#' \enumerate{
#'   \item \code{normalize01}: Normalize continuous attributes to [0,1]  
#'   \item \code{Cons.Vars}: Detect binary variables with extreme imbalance in 
#' class distribution 
#'   \item \code{opt_thresh}: Wrapper for the function to compute optimal classification threshold 
#' from the package \code{\link[PresenceAbsence]{PresenceAbsence}} 
#'   \item \code{Performance.measures}: compute performance measures: AUC, 
#' sensitivity, specificity, PPV, etc.
#'   \item \code{get.para}:  
#'   Extract tuning parameters from cross-validation (RunCV) model with best performance 
#'   \item \code{get.cvResults}: Extract cross-validation results: average
#'    performance measures over all cross-validations 
#' \item \code{get.bootResults}: Extract bootstrap results: mean and confidence 
#' intervals 
#' \item \code{get.varimp}: Extract variable importance for GLM, GBM and RF models  
#' }
#'
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#'
#' @importFrom PresenceAbsence optimal.thresholds
#' @importFrom PresenceAbsence presence.absence.accuracy
#' @importFrom PresenceAbsence predicted.prevalence
#' @importFrom plyr ddply numcolwise dlply
#'
#'
NULL 
#'
#' @rdname  Utils
#' @export
Split <- function(Nt, K){
size = Nt%/%K
A = runif(Nt)
rk = rank(A)
bk = factor((rk-1)%/%size + 1)
if(Nt%%K > 0) {
levels(bk)[levels(bk) == K+1]  = 1 
}
return(bk)
}
#' @rdname  Utils
#' @export
SplitCrossVal <- function(D.dat, k, ix.cv, val.propo = 0.15){
dat.tmp <- D.dat[k!= ix.cv, ,drop=FALSE]
dat.tst <- D.dat[k==ix.cv, ,drop=FALSE]
ix <- sample(nrow(dat.tmp), floor(nrow(dat.tmp)*val.propo))
dat.val <-  dat.tmp[ix, ,drop=FALSE] 
dat.trn <-  dat.tmp[-ix, ,drop=FALSE] 
res <- list(dat.trn = dat.trn, dat.val = dat.val, dat.tst = dat.tst)
return(res)
}
#' @rdname  Utils
#' @export
collect.garbage = function(){
#if (exists("collect.garbage") ) rm(collect.garbage)
## The following function collects garbage until the memory is clean.
## Usage: 1. immediately call this function after you call a function or
##        2. rm()
	while (gc()[2,4] != gc()[2,4]){}
}
#' @rdname  Utils
#' @export
CreateDummy = function(D.dat) {
nme <- colnames(D.dat) 
res <- do.call(cbind.data.frame, lapply(nme, function(x) {
  A <- as.character(D.dat[,x])
  levs <- sort(unique(A))
  Cols <- match(A, levs)
  m <- matrix(0, nrow = length(D.dat[,x]), ncol = length(levs),
              dimnames = list(NULL, paste(x, levs, sep = ".")))
  m[cbind(sequence(length(D.dat[,x])), Cols)] <- 1L
m
}))
return(res)
}

#' @rdname  Utils
#' @export
# normalize data to (0,1)
normalize01 = function(x){
doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - 
min(x, na.rm=TRUE))} 
as.data.frame(lapply(x, doit))  
}
#'
#' @rdname  Utils
#' @export
Cons.Vars = function(dummy, thresh=0.00016, uniqueCut = 10){
lunique = length(unique(dummy))
nme <- FALSE 
if(lunique == 1){
nme <- TRUE
}
if(lunique == 2){
levs <- sort(unique(dummy))
m1 = mean(dummy ==levs[1], na.rm = TRUE)
m0 = mean(dummy ==levs[2], na.rm = TRUE)
if(m1 <= thresh | m0 <= thresh) nme <- TRUE
} 
return(nme)
}
#' @rdname  Utils
#' @export 
opt_thresh <- function(prob, obs){
thresh = 0.5 
if(length(unique(obs)) > 1){
obs <- as.numeric(as.factor(obs))-1 
SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = prob)
thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = 9)
thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
}
return(thresh)
}
#' @rdname  Utils 
#' @export
Performance.measures <- function(pred, obs, prevalence = NULL, threshold=NULL){
obs <- as.numeric(factor(obs))-1 
## get best cut-off 
if(is.null(threshold))
threshold <- opt_thresh(prob=pred, obs=obs)
### get these performance measures
nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
"specificity", "specificity.sd")
xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
pred.prev <- predicted.prevalence(DATA=xx, threshold = threshold)[, 
c("Obs.Prevalence", "Predicted")]
nme <- c("Pos Pred Value", "Neg Pred Value", "Balanced Accuracy")
accuracy$G.mean <- sqrt(as.numeric(accuracy$sensitivity)*as.numeric(accuracy$specificity))
accuracy$BER <- 1 - 0.5*(as.numeric(accuracy$sensitivity) + as.numeric(accuracy$specificity))  
if(is.null(prevalence)) prevalence = as.numeric(pred.prev$Obs.Prevalence)
obs <- factor(ifelse(obs == 1, "Yes", "No"), levels = c("Yes", "No"))
pred <- factor(ifelse(pred >= threshold, "Yes", "No"), levels = c("Yes", "No"))
cmx <- confusionMatrix(data=pred, reference=obs,  prevalence = prevalence)$byClass[nme]
res <- cbind.data.frame(accuracy, t(cmx), pred.prev, threshold = threshold)
return(res)
}
#'
#' @rdname  Utils
#' @export
get.para <- function(object, para){
resp.vars <- object$resp.vars 
pp <- lapply(resp.vars, function(zz){ 

tab <- do.call(rbind.data.frame, lapply(object$CV.res, function(x) {   ## each CV 
cls <- names(x[[zz]]) 
do.call(rbind.data.frame, lapply(1:length(x[[zz]]), function(y){    ## each classifier 
cbind.data.frame(classifier = cls[y], AUC = x[[zz]][[y]]$perf$tst.perf$AUC)
}))}))

ix.tab = ddply(tab, .variables = "classifier", .fun = numcolwise(which.max))
ix <- ix.tab$AUC 
cls <- ix.tab$classifier
for(ii in 1:length(ix))
para[[cls[ii]]] <- object$CV.res[[ix[ii]]][[zz]][[cls[ii]]]$para[[cls[ii]]]
para
}
)
names(pp) <- resp.vars
pp
}
#'
#' @rdname  Utils
#' @export
get.cvResults <- function(object){
tab <- do.call(rbind.data.frame, lapply(object$CV.res, function(x) {   ## each CV 
models <- names(x) 
do.call(rbind.data.frame, lapply(models,  function(y){
cls <- names(x[[y]])
do.call(rbind.data.frame, lapply(cls, function(z){
cbind.data.frame(Outcome=y, Classifier = z, x[[y]][[z]]$perf$tst.perf)
}))}))}))

ddply(tab, .variables = c("Outcome","Classifier"), .fun = numcolwise(mean) )[, -c(3,4)]
}
#'
#' @rdname  Utils
#' @export
get.bootResults <- function(object, alpha=0.05){

tab <- do.call(rbind.data.frame, lapply(object$Boot.res, function(x) {   ## each CV 
models <- names(x) 
do.call(rbind.data.frame, lapply(models,  function(y){
cls <- names(x[[y]])
do.call(rbind.data.frame, lapply(cls, function(z){
cbind.data.frame(Outcome=y, Classifier = z, x[[y]][[z]]$perf$tst.perf)
}))}))}))

MNE <- ddply(tab, .variables = c("Outcome","Classifier"), .fun = numcolwise(mean))[, -c(3,4)]
CI <-   ddply(tab, .variables = c("Outcome","Classifier"),  numcolwise(quantile), 
                  probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)[, -c(3,4)]
list(Mean = MNE, CI = CI)                
}

#'
#' @rdname  Utils
#' @export
get.varimp <- function(object, alpha = 0.05, top = 10){
if(!object$varimp) stop("no variable importance available in model") 

cls <- c("GLM", "GBM", "RF")
tab <- do.call(rbind.na, lapply(object$CV.res, function(x) {   ## each CV 
models <- names(x) 
do.call(rbind.na, lapply(models,  function(y){  ## outcome 

GLM.vars = rownames(x[[y]][["GLM"]]$varimp); GLM.imp = x[[y]][["GLM"]]$varimp[, "importance"] 
GBM.vars = rownames(x[[y]][["GBM"]]$varimp); GBM.imp = x[[y]][["GBM"]]$varimp[, "importance"] 
A <- min(GBM.imp); B <- max(GBM.imp) ; 	a = 0; b = 100; GBM.imp = (((GBM.imp - A)*(b-a))/(B-A) + a)        
RF.vars = rownames(x[[y]][["RF"]]$varimp); RF.imp = x[[y]][["RF"]]$varimp[, "importance"]
A <- min(RF.imp); B <- max(RF.imp) ; 	a = 0; b = 100; RF.imp = (((RF.imp - A)*(b-a))/(B-A) + a)        
cbind.na(Outcome = rep(y, max(length(GLM.vars), length(GBM.vars),length(RF.vars) )), 
        GLM.vars = GLM.vars, GLM.imp = GLM.imp, GBM.vars = GBM.vars,
         GBM.imp =GBM.imp, RF.vars = RF.vars, RF.imp = RF.imp)
} ))}))
tab <- data.frame(tab)
tab <- within(tab, {
GLM.imp = as.numeric(GLM.imp)
GBM.imp = as.numeric(GBM.imp)
RF.imp = as.numeric(RF.imp)})


vimp = dlply(tab, .variables = "Outcome", .fun = function(x) x)

XX.vimp <- lapply(vimp, function(x){
xx <- ddply(x[, c("GLM.vars", "GLM.imp")], .variables = c("GLM.vars"), numcolwise(mean), na.rm = TRUE) 
GLM <-  xx[order(xx$GLM.imp, decreasing = TRUE), ,drop=FALSE]     
rownames(GLM) <- NULL                        
xx <- ddply(x[, c("GBM.vars", "GBM.imp")], .variables = c("GBM.vars"), numcolwise(mean), na.rm = TRUE)
GBM <- xx[order(xx$GBM.imp, decreasing = TRUE), ,drop=FALSE]    
rownames(GBM) <- NULL                        
xx <- ddply(x[, c("RF.vars", "RF.imp")], .variables = c("RF.vars"), numcolwise(mean), na.rm = TRUE) 
RF <- xx[order(xx$RF.imp, decreasing = TRUE), ,drop=FALSE]     
rownames(RF) <- NULL 
            
GLM[,-c(1)]  <- format(round(GLM[, -c(1)], 2), nsmall=2)
GBM[,-c(1)]  <- format(round(GBM[, -c(1)], 2), nsmall=2)
RF[,-c(1)]  <- format(round(RF[, -c(1)], 2), nsmall=2)
tb <- cbind.na(GLM, GBM[, ,drop=FALSE], RF[, ,drop=FALSE])
if(nrow(tb) > top) tb <- tb[1:top, ] 
rownames(tb) <- NULL
tb
})
XX.vimp
}

#' @rdname  Utils 
#' @export
cbind.na <- function (..., deparse.level = 1)
{
 na <- nargs() - (!missing(deparse.level))
 deparse.level <- as.integer(deparse.level)
 stopifnot(0 <= deparse.level, deparse.level <= 2)
 argl <- list(...)
 while (na > 0 && is.null(argl[[na]])) {
 argl <- argl[-na]
 na <- na - 1
 }
 if (na == 0)
 return(NULL)
 if (na == 1) {
 if (isS4(..1))
 return(cbind2(..1))
 else return(matrix(...)) ##.Internal(cbind(deparse.level, ...)))
 }
 if (deparse.level) {
 symarg <- as.list(sys.call()[-1L])[1L:na]
 Nms <- function(i) {
 if (is.null(r <- names(symarg[i])) || r == "") {
 if (is.symbol(r <- symarg[[i]]) || deparse.level ==
 2)
 deparse(r)
 }
 else r
 }
 }
 ## deactivated, otherwise no fill in with two arguments
 if (na == 0) {
 r <- argl[[2]]
 fix.na <- FALSE
 }
 else {
 nrs <- unname(lapply(argl, nrow))
 iV <- sapply(nrs, is.null)
 fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
 ## deactivated, otherwise data will be recycled
 #if (fix.na) {
 # nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
 # argl[[na]] <- cbind(rep(argl[[na]], length.out = nr),
 # deparse.level = 0)
 #}
 if (deparse.level) {
 if (fix.na)
 fix.na <- !is.null(Nna <- Nms(na))
 if (!is.null(nmi <- names(argl)))
 iV <- iV & (nmi == "")
 ii <- if (fix.na)
 2:(na - 1)
 else 2:na
 if (any(iV[ii])) {
 for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
 names(argl)[i] <- nmi
 }
 }

 ## filling with NA's to maximum occuring nrows
 nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
 maxRow <- max(nRow, na.rm = TRUE)
 argl <- lapply(argl, function(x) if (is.null(nrow(x))) c(x, rep(NA, maxRow - length(x)))
 else rbind.na(x, matrix(, maxRow - nrow(x), ncol(x))))
 r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
 }
 d2 <- dim(r)
 r <- cbind2(argl[[1]], r)
 if (deparse.level == 0)
 return(r)
 ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
 ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
 if (ism1 && ism2)
 return(r)
 Ncol <- function(x) {
 d <- dim(x)
 if (length(d) == 2L)
 d[2L]
 else as.integer(length(x) > 0L)
 }
 nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
 nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
 if (nn1 || nn2 || fix.na) {
 if (is.null(colnames(r)))
 colnames(r) <- rep.int("", ncol(r))
 setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams))
 ""
 else nams
 if (nn1)
 setN(1, N1)
 if (nn2)
 setN(1 + l1, N2)
 if (fix.na)
 setN(ncol(r), Nna)
 }
 r
}
 
#'
#' @rdname  Utils
#' @export
rbind.na <- function (..., deparse.level = 1)
{
    na <- nargs() - (!missing(deparse.level))
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)
    argl <- list(...)
    while (na > 0 && is.null(argl[[na]])) {
        argl <- argl[-na]
        na <- na - 1
    }
    if (na == 0)
        return(NULL)
    if (na == 1) {
        if (isS4(..1))
            return(rbind2(..1))
        else return(matrix(..., nrow = 1)) ##.Internal(rbind(deparse.level, ...)))
    }
    if (deparse.level) {
        symarg <- as.list(sys.call()[-1L])[1L:na]
        Nms <- function(i) {
            if (is.null(r <- names(symarg[i])) || r == "") {
                if (is.symbol(r <- symarg[[i]]) || deparse.level ==
                  2)
                  deparse(r)
            }
            else r
        }
    }

    ## deactivated, otherwise no fill in with two arguments
    if (na == 0) {
        r <- argl[[2]]
        fix.na <- FALSE
    }
    else {
        nrs <- unname(lapply(argl, ncol))
        iV <- sapply(nrs, is.null)
        fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
        ## deactivated, otherwise data will be recycled
        #if (fix.na) {
        #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
        #    argl[[na]] <- rbind(rep(argl[[na]], length.out = nr),
        #        deparse.level = 0)
        #}
        if (deparse.level) {
            if (fix.na)
                fix.na <- !is.null(Nna <- Nms(na))
            if (!is.null(nmi <- names(argl)))
                iV <- iV & (nmi == "")
            ii <- if (fix.na)
                2:(na - 1)
            else 2:na
            if (any(iV[ii])) {
                for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i)))
                  names(argl)[i] <- nmi
            }
        }

        ## filling with NA's to maximum occuring ncols
        nCol <- as.numeric(sapply(argl, function(x) if (is.null(ncol(x))) length(x)
                                                    else ncol(x)))
        maxCol <- max(nCol, na.rm = TRUE)
        argl <- lapply(argl, function(x)  if (is.null(ncol(x))) c(x, rep(NA, maxCol - length(x)))
                                          else cbind(x, matrix(, nrow(x), maxCol - ncol(x))))

        ## create a common name vector from the
        ## column names of all 'argl' items
        namesVEC <- rep(NA, maxCol)
        for (i in 1:length(argl)) {
          CN <- colnames(argl[[i]])
          m <- !(CN %in% namesVEC)
          namesVEC[m] <- CN[m]
        }

        ## make all column names from common 'namesVEC'
        for (j in 1:length(argl)) {
          if (!is.null(ncol(argl[[j]]))) colnames(argl[[j]]) <- namesVEC
        }

        r <- do.call(rbind, c(argl[-1L], list(deparse.level = deparse.level)))
    }

    d2 <- dim(r)

    ## make all column names from common 'namesVEC'
    colnames(r) <- colnames(argl[[1]])

    r <- rbind2(argl[[1]], r)

    if (deparse.level == 0)
        return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
    if (ism1 && ism2)
        return(r)
    Nrow <- function(x) {
        d <- dim(x)
        if (length(d) == 2L)
            d[1L]
        else as.integer(length(x) > 0L)
    }
    nn1 <- !is.null(N1 <- if ((l1 <- Nrow(..1)) && !ism1) Nms(1))
    nn2 <- !is.null(N2 <- if (na == 2 && Nrow(..2) && !ism2) Nms(2))
    if (nn1 || nn2 || fix.na) {
        if (is.null(rownames(r)))
            rownames(r) <- rep.int("", nrow(r))
        setN <- function(i, nams) rownames(r)[i] <<- if (is.null(nams))
            ""
        else nams
        if (nn1)
            setN(1, N1)
        if (nn2)
            setN(1 + l1, N2)
        if (fix.na)
            setN(nrow(r), Nna)
    }
    r
}
 
#'
#' @rdname  Utils
#' @export
data.frame.na <- function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,
    stringsAsFactors = FALSE)
{
    data.row.names <- if (check.rows && is.null(row.names))
        function(current, new, i) {
            if (is.character(current))
                new <- as.character(new)
            if (is.character(new))
                current <- as.character(current)
            if (anyDuplicated(new))
                return(current)
            if (is.null(current))
                return(new)
            if (all(current == new) || all(current == ""))
                return(new)
            stop(gettextf("mismatch of row names in arguments of 'data.frame', item %d",
                i), domain = NA)
        }
    else function(current, new, i) {
        if (is.null(current)) {
            if (anyDuplicated(new)) {
                warning("some row.names duplicated: ", paste(which(duplicated(new)),
                  collapse = ","), " --> row.names NOT used")
                current
            }
            else new
        }
        else current
    }
    object <- as.list(substitute(list(...)))[-1L]
    mrn <- is.null(row.names)
    x <- list(...)
    n <- length(x)
    if (n < 1L) {
        if (!mrn) {
            if (is.object(row.names) || !is.integer(row.names))
                row.names <- as.character(row.names)
            if (any(is.na(row.names)))
                stop("row names contain missing values")
            if (anyDuplicated(row.names))
                stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                  collapse = ", "))
        }
        else row.names <- integer(0L)
        return(structure(list(), names = character(0L), row.names = row.names,
            class = "data.frame"))
    }
    vnames <- names(x)
    if (length(vnames) != n)
        vnames <- character(n)
    no.vn <- !nzchar(vnames)
    vlist <- vnames <- as.list(vnames)
    nrows <- ncols <- integer(n)
    for (i in seq_len(n)) {
        xi <- if (is.character(x[[i]]) || is.list(x[[i]]))
            as.data.frame(x[[i]], optional = TRUE, stringsAsFactors = stringsAsFactors)
        else as.data.frame(x[[i]], optional = TRUE)
        nrows[i] <- .row_names_info(xi)
        ncols[i] <- length(xi)
        namesi <- names(xi)
        if (ncols[i] > 1L) {
            if (length(namesi) == 0L)
                namesi <- seq_len(ncols[i])
            if (no.vn[i])
                vnames[[i]] <- namesi
            else vnames[[i]] <- paste(vnames[[i]], namesi, sep = ".")
        }
        else {
            if (length(namesi))
                vnames[[i]] <- namesi
            else if (no.vn[[i]]) {
                tmpname <- deparse(object[[i]])[1L]
                if (substr(tmpname, 1L, 2L) == "I(") {
                  ntmpn <- nchar(tmpname, "c")
                  if (substr(tmpname, ntmpn, ntmpn) == ")")
                    tmpname <- substr(tmpname, 3L, ntmpn - 1L)
                }
                vnames[[i]] <- tmpname
            }
        }
        if (missing(row.names) && nrows[i] > 0L) {
            rowsi <- attr(xi, "row.names")
            nc <- nchar(rowsi, allowNA = FALSE)
            nc <- nc[!is.na(nc)]
            if (length(nc) && any(nc))
                row.names <- data.row.names(row.names, rowsi,
                  i)
        }
        nrows[i] <- abs(nrows[i])
        vlist[[i]] <- xi
    }
    nr <- max(nrows)
    for (i in seq_len(n)[nrows < nr]) {
        xi <- vlist[[i]]
        if (nrows[i] > 0L) {
            xi <- unclass(xi)
            fixed <- TRUE
            for (j in seq_along(xi)) {
                ### added NA fill to max length/nrow
                xi1 <- xi[[j]]
                if (is.vector(xi1) || is.factor(xi1))
                  xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                else if (is.character(xi1) && class(xi1) == "AsIs")
                  xi[[j]] <- structure(c(xi1, rep(NA, nr - nrows[i])),
                    class = class(xi1))
                else if (inherits(xi1, "Date") || inherits(xi1,
                  "POSIXct"))
                  xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                else {
                  fixed <- FALSE
                  break
                }
            }
            if (fixed) {
                vlist[[i]] <- xi
                next
            }
        }
        stop("arguments imply differing number of rows: ", paste(unique(nrows),
            collapse = ", "))
    }
    value <- unlist(vlist, recursive = FALSE, use.names = FALSE)
    vnames <- unlist(vnames[ncols > 0L])
    noname <- !nzchar(vnames)
    if (any(noname))
        vnames[noname] <- paste("Var", seq_along(vnames), sep = ".")[noname]
    if (check.names)
        vnames <- make.names(vnames, unique = TRUE)
    names(value) <- vnames
    if (!mrn) {
        if (length(row.names) == 1L && nr != 1L) {
            if (is.character(row.names))
                row.names <- match(row.names, vnames, 0L)
            if (length(row.names) != 1L || row.names < 1L ||
                row.names > length(vnames))
                stop("row.names should specify one of the variables")
            i <- row.names
            row.names <- value[[i]]
            value <- value[-i]
        }
        else if (!is.null(row.names) && length(row.names) !=
            nr)
            stop("row names supplied are of the wrong length")
    }
    else if (!is.null(row.names) && length(row.names) != nr) {
        warning("row names were found from a short variable and have been discarded")
        row.names <- NULL
    }
    if (is.null(row.names))
        row.names <- .set_row_names(nr)
    else {
        if (is.object(row.names) || !is.integer(row.names))
            row.names <- as.character(row.names)
        if (any(is.na(row.names)))
            stop("row names contain missing values")
        if (anyDuplicated(row.names))
            stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                collapse = ", "))
    }
    attr(value, "row.names") <- row.names
    attr(value, "class") <- "data.frame"
    value
}




