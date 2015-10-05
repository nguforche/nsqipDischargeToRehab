#' Train, Validate and Test various models  
#' 
#' Functions to train, select best hyper-paramters, select best classification thresholds 
#'  and evaluate final model on the test set  using a number of machine 
#' learning classification algorithms 
#' @name Train 
#' @param  classifier character list of classification models. See names(TrainAllModels()). 
#' @param  X.trn,X.val,X.tst matrix of predictors 
#' @param  Y.trn,Y.val,Y.tst matrix of binary \{0,1\} response variable
#' @param  varimp (logical) compute variable importance ?
#' @param  opt.para (logical) tune parameters ?
#' @param  para named list of model parameters 
#' @param  return.model (logical) return trained model ?  
#' @param \dots further arguments passed to or from other methods.
#' @return \code{TrainModels} returns a list of functions for training various algorithms. Currently 
#' 7 machine learning algorithms ELR, GLM, GLMnet, RF, GBM, avNNET, and SVM. Others can be easily added. 
#' \code{Train.Validate.Test} returns a list of performance measures for each classifier traoned:
#' \item{model}{trained model}
#' \item{para}{named list of model hyper-paramters (tunned values if opt.para = TRUE)}
#' \item{run.time}{compute time} 
#' \item{varimp}{variable importance if classifier is GLM, GBM or RF}
#' \item{perf}{a list with two data.frames: val.perf and tst.perf containing performance measures 
#' for validation and test sets }

#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#' @import ELR caret kernlab 
#'
NULL 
#' @rdname Train 
#' @export
#' @examples
#' Trainmodels <- lapply(TrainModels(), function(x) x) 
#'
TrainModels <- function(){
res = list(
ELR = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){ 
	result <- vector("list", length = 5)
	names(result) <- c("trn", "val", "tst", "main.model", "para")
	class(result) <- "ELR"	
	result$para <- para 		
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 		
    dt <- cbind.data.frame(resp= Y[,1], X)
    colnames(dt) <- c(resp, rhs.vars) 
    form <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))
         
		if(opt.para) {  
		grid <- para$ELR$grid 
		g.mean <- mod1 <- c() 
		
		for(ii in 1:nrow(grid)){
		    para$p <- grid[ii, "p"]
		    para$gamma <- grid[ii, "gamma"]         
		    mod <-  ELR(form, dt, para$ELR)        
		    pp <- predict(mod, X.val)[,2]		
		    pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))
		    g.mean <- c(g.mean, as.numeric(Performance.measures(pred=pp, obs=X.val[, resp])[, "G.mean"]))                
		    mod1[[ii]] <- mod 
	    }	
	           ix = which.max(g.mean)
		  if(is.na(ix)){
		   result$main.model <- ELR(form, dt,  para)
		  } else {
			result$main.model <- mod1[[ix]]
			result$para$ELR$p <- grid[ix, "p"]
			result$para$ELR$gamma <- grid[ix, "gamma"] 
		   }
	  } else {
	    result$main.model <- ELR(form, dt,  para$ELR)
		}
			
		pp <- predict(result$main.model, X)[, 2]	
		pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
		result$trn <- pp		
		pp <- predict(result$main.model, X.val[, rhs.vars])[, 2]
		pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
		result$val <- pp		
		pp <- predict(result$main.model, X.tst)[, 2]
		pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
		result$tst <- pp 
		result$form <- form 
    return(result)
},

GLM = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <-c("trn", "val", "tst", "main.model", "para") 
	class(result) <- "GLM"
	result$para <- para 
			
	rhs.vars <- colnames(X)
	resp <- colnames(Y)	

	fitControl <- trainControl(method = "none", classProbs = TRUE)
	
	dep.vars <- findLinearCombos(X)$remove
	
	if(length(dep.vars) < length(rhs.vars)-1) 
	rhs.vars <- rhs.vars[-dep.vars]
	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      
	
	dt <- cbind.data.frame(resp= Y[,1], X[, rhs.vars, drop=FALSE])
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))
	colnames(dt) <- c(resp, rhs.vars)
	
	result$main.model <-  train(formula.string, data = dt, method = "glm", family = "binomial", 
                                     trControl = fitControl, metric = "ROC", control = list(maxit = 100))     
	pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[,2]        
	result$trn <- pp          
	pp <-  predict(result$main.model, newdata = X.val[, rhs.vars], type = "prob")[,2]        
	result$val <- pp                                                           
	pp <-  predict(result$main.model, newdata = X.tst[, rhs.vars], type = "prob")[,2] 
	result$tst <- pp 
	result$form <- formula.string  
	return(result)
},

GLMnet = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <-c("trn", "val", "tst", "main.model", "para") 
	class(result) <- "GLMnet"
	result$para <- para 			
	rhs.vars <- colnames(X)
	resp <- colnames(Y)	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))
	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))
	colnames(dt) <- c(resp, rhs.vars)	
		      
	if(opt.para) {
	     grid <- para$GLMnet$grid	
	     g.mean <- mod1 <- c() 
		 for(ii in 1:nrow(grid)){ 
			fitControl <- trainControl(method = "none")	
			mod <-  train(formula.string, data = dt, method = "glmnet", 
				                  family = "binomial", 
				                  tuneGrid = data.frame(alpha = grid[ii, "alpha"], 
				                  lambda = grid[ii, "lambda"]), trControl = fitControl)		 
                                                                             
			pp <- predict(mod, newdata = X.val[, rhs.vars], type = "prob")[,2]
			pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
		        g.mean <- c(g.mean, as.numeric(Performance.measures(pred=pp, obs=X.val[, resp])[, "G.mean"]))                
			mod1[[ii]] <- mod 
	     }	
	    ix = which.max(g.mean)
        if(is.na(ix)){
		    fitControl <- trainControl(method = "none")	
		    result$main.model <-  train(formula.string, data = dt, method = "glmnet", 
		                              family = "binomial", tuneGrid = data.frame(alpha = para$GLMnet$alpha, 
		                              lambda = para$GLMnet$lambda), trControl = fitControl)
        } else {		     
		result$main.model <- mod1[[ix]]
		result$para$GLMnet$alpha  <-  grid[ix, "alpha"]
		result$para$GBMnet$lambda  <-  grid[ix, "lambda"]
        }
 	} else {
		    fitControl <- trainControl(method = "none")	
		    result$main.model <-  train(formula.string, data = dt, method = "glmnet", 
		                              family = "binomial", tuneGrid = data.frame(alpha = para$GLMnet$alpha, 
		                              lambda = para$GLMnet$lambda), trControl = fitControl)
        }                                                 
        pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[,2]        
        result$trn <- pp          
        pp <-  predict(result$main.model, newdata = X.val[, rhs.vars], type = "prob")[,2]        
        result$val <- pp                                                           
        pp <-  predict(result$main.model, newdata = X.tst[, rhs.vars], type = "prob")[,2] 
        result$tst <- pp   
	    result$form <- formula.string           
	return(result)
},
#### random forest  
RF = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <- c("trn", "val", "tst", "main.model", "para")
	class(result) <- "RF"	
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 
	result$para <- para 	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      
	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))	
	colnames(dt) <- c(resp, rhs.vars)
		
	if(opt.para) { 
         mtry <- floor(seq(1, floor(length(rhs.vars)/2), 2.5)) 
	     fitControl <- trainControl(method = "none", verboseIter = FALSE)
	     g.mean <- mod1 <- c() 
		 for(ii in 1:length(mtry)){ 
			 mod <- train(X, dt[, resp],  method = "rf", trControl = fitControl,
			                      tuneGrid = data.frame(mtry = mtry[ii] ), 
			                      ntree = para$RF$ntree, importance = FALSE, proximity=FALSE)	
			                                                                         
		    pp <- predict(mod, newdata = X.val[, rhs.vars], type = "prob")[,2]
		    pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
		    g.mean <- c(g.mean, as.numeric(Performance.measures(pred=pp, obs=X.val[, resp])[, "G.mean"]))                
	            mod1[[ii]] <- mod 
	     }	
	    ix = which.max(g.mean)
	    if(is.na(ix)){
			fitControl <- trainControl(method = "none", verboseIter = FALSE)	        
			result$main.model <- train(X, dt[, resp],  method = "rf", trControl = fitControl,
			                      tuneGrid = data.frame(mtry = floor(ncol(X)/3)), 
			                      ntree = para$RF$ntree, importance = FALSE, proximity=FALSE)		      
	      } else {	     
			result$main.model <- mod1[[ix]]
			result$para$RF$mtry <- mtry[ix]
	    }
	
  } else{
	    fitControl <- trainControl(method = "none", verboseIter = FALSE)	        
	    result$main.model <- train(X, dt[, resp],  method = "rf", trControl = fitControl,
	                          tuneGrid = data.frame(mtry = para$RF$mtry), 
	                          ntree = para$RF$ntree, importance = FALSE, proximity=FALSE)	
         }    
        pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
        result$trn <- pp                       
        pp <-  predict(result$main.model, newdata = X.val[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
        result$val <- pp
        pp <- predict(result$main.model, newdata =X.tst[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
        result$tst <- pp  
        result$form <- formula.string
	return(result)
},
#### GBM   
GBM = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <- c("trn", "val", "tst", "main.model", "para")
	class(result) <- "GBM"	
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 
	result$para <- para 	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      
	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))	
	colnames(dt) <- c(resp, rhs.vars)
		
	if(opt.para) { 
	     grid <- para$GBM$grid 
	     fitControl <- trainControl(method = "none", classProbs = TRUE)
	     g.mean <- mod1 <- c() 
		 for(ii in 1:nrow(grid)){ 
		     mod <-  train(formula.string, dt, method = "gbm", trControl = fitControl,
		                       verbose = FALSE, tuneGrid = data.frame(n.trees = grid[ii, "n.trees"], 
		                       interaction.depth=grid[ii, "interaction.depth"], 
		                       shrinkage=grid[ii, "shrinkage"], n.minobsinnode=grid[ii, "n.minobsinnode"]), 
		                       metric = "ROC")
		                                                                     
		     pp <-  predict(mod, newdata = X.val, type = "prob")[,2] 
		     pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
		     g.mean <- c(g.mean, as.numeric(Performance.measures(pred=pp, obs=X.val[, resp])[, "G.mean"]))                
		     mod1[[ii]] <- mod 
	    }
	     ix = which.max(g.mean)
	     if(is.na(ix)){
			 fitControl <- trainControl(method = "none", classProbs = TRUE)
			 result$main.model <-  train(formula.string, dt, method = "gbm", trControl = fitControl,
		                                verbose = FALSE, 
		                                tuneGrid = data.frame(n.trees = para$GBM$n.trees, 
		                                interaction.depth=para$GBM$interaction.depth, 
		                                shrinkage=para$GBM$shrinkage, 
		                                n.minobsinnode=para$GBM$n.minobsinnode), 
		                                metric = "ROC")	      
	      } else {	
		    result$main.model <- mod1[[ix]]
		    result$para$GBM$n.trees <- grid[ix, "n.trees"]
		    result$para$GBM$interaction.depth  <-  grid[ix, "interaction.depth"]
		    result$para$GBM$shrinkage  <- grid[ix, "shrinkage"]
		    result$para$GBM$n.minobsinnode  <- grid[ix, "n.minobsinnode"]
        }                                    
   } else{
	     fitControl <- trainControl(method = "none", classProbs = TRUE)
	     result$main.model <-  train(formula.string, dt, method = "gbm", trControl = fitControl,
                                    verbose = FALSE, 
                                    tuneGrid = data.frame(n.trees = para$GBM$n.trees, 
                                    interaction.depth=para$GBM$interaction.depth, 
                                    shrinkage=para$GBM$shrinkage, 
                                    n.minobsinnode=para$GBM$n.minobsinnode), 
                                    metric = "ROC")
         } 
        pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
        result$trn <- pp                          
        pp <-  predict(result$main.model, newdata = X.val[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))                 
        result$val <- pp
        pp <- predict(result$main.model, newdata =X.tst[, rhs.vars], type = "prob")[,2] 
        pp <- ifelse( abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))  
        result$tst <- pp  
        result$form <- formula.string
	return(result)
},
## penalized LDA 
pLDA = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <- c("trn", "val", "tst", "main.model", "para")
	class(result) <- "pFDA"	
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 
	result$para <- para 	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      
	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))	
	colnames(dt) <- c(resp, rhs.vars)
		
    if(opt.para) { 
      fitControl <- trainControl(method = para$method, number =  para$number)                           	                                 	
      result$main.model <-  train(X, dt[, resp],  method = "pda", trControl = fitControl,
                                  tuneLength = para$tuneLength)       
      result$para$pLDA$lambda  <-  result$main.model$bestTune$lambda
    } else{
      fitControl <- trainControl(method = "none")
      result$main.model <-  train(X, dt[, resp],  method = "pda", trControl = fitControl,
                                  tuneGrid = data.frame(lambda = para$pLDA$lambda))
    } 
    pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[,2] 
    result$trn <- pp                          
    pp <-  predict(result$main.model, newdata = X.val[, rhs.vars], type = "prob")[,2] 
    result$val <- pp
    pp <- predict(result$main.model, newdata =X.tst[, rhs.vars], type = "prob")[,2] 
    result$tst <- pp  
    return(result)
},
#### average neural networks  
avNNET = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <- c("trn", "val", "tst", "main.model", "para")
	class(result) <- "avNNET"	
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 
	result$para <- para 	
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      

	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))	
	colnames(dt) <- c(resp, rhs.vars)	
	if(opt.para) { 
	fitControl <- trainControl(method = para$method, number = para$number, verboseIter = FALSE)	
	result$main.model <-  train(X, dt[, resp],  method = "avNNet", trControl = fitControl,
	                        tuneLength = para$tuneLength, maxit = para$avNNET$maxit,  
                                    MaxNWts = 10000, trace = FALSE)                          
                           
        result$para$avNNET$nnet.size  <-  result$main.model$bestTune$size 
        result$para$avNNET$nnet.decay <-  result$main.model$bestTune$decay 
        result$para$avNNET$nnet.bag   <-  result$main.model$bestTune$bag               
        } else{
	fitControl <- trainControl(method = "none")
	result$main.model <-  train(X, dt[, resp],  method = "nnet", 
	                       trControl = fitControl,
                               verbose = FALSE, tuneGrid = 
                               data.frame(size = para$avNNET$nnet.size, 
                               decay=para$avNNET$nnet.decay), 
                               maxit = para$avNNET$maxit, MaxNWts = 10000, 
                               trace = FALSE)                          	
	}                   
        pp <-  predict(result$main.model, newdata = X[, rhs.vars], type = "prob")[, 2] 
        result$trn <- pp
        pp <- predict(result$main.model, newdata =X.val[, rhs.vars], type = "prob")[,2] 
        result$val <- pp  
        pp <- predict(result$main.model, newdata =X.tst[, rhs.vars], type = "prob")[,2] 
        result$tst <- pp          
	return(result)
},
#### support vector machines radial 
SVM = function(X, Y, X.val, X.tst, para, opt.para = FALSE, ...){
	result <- vector("list", length = 5)
	names(result) <-  c("trn", "val", "tst", "main.model", "para")
	class(result) <- "SVM"	
	rhs.vars <- colnames(X)
	resp <- colnames(Y) 
	result$para <- para 
	formula.string <- as.formula(paste0(paste0(resp, " ~"), paste0(rhs.vars, collapse= "+")))      

	dt <- cbind.data.frame(resp= Y[,1], X)
	dt[, "resp"] <- factor(ifelse(dt[, "resp"] == 1, "Yes", "No"))	
	colnames(dt) <- c(resp, rhs.vars)
	
	if(opt.para) { 
	    grid <- para$SVM$grid
        n = nrow(grid) 
        AUC = c()
        mod <- lapply(rep("ksvm", n), new )
        for(ii in 1:n){
	  	    mod[[ii]] <-  ksvm(formula.string, data = dt,  type = "C-bsvc", 
	  	                  kernel = para$SVM$kernel,  C =  grid[ii, "C"], 
	  	                  kpar = list(sigma = grid[ii, "sigma"]), cross = 0, 
	  	                  prob.model = TRUE, scaled = TRUE)
		           
		     pp = abs(predict(mod[[ii]], newdata=X, type="probabilities"))[,2]
		    AUC <- c(AUC, as.numeric(Performance(pp, Y[,1])$AUC))
        }
		ix = which.max(AUC)
		if(is.na(ix)){	
	                  result$main.model <- ksvm(formula.string, data = dt,  type = "C-bsvc", 
	                                       kernel = para$SVM$kernel, C =  para$SVM$C, 
	                                       kpar = list(sigma = para$SVM$sigma), cross = 0, 
	                                       prob.model = TRUE, scaled = TRUE)				
		} else {
		result$main.model <- mod[[ix]]
		result$para$SVM$C <-  grid[ix, "C"]
		result$para$SVM$sigma <-  grid[ix, "sigma"]
		}		
   } else{
	result$main.model <- ksvm(formula.string, data = dt,  type = "C-bsvc", 
	                     kernel = para$SVM$kernel, C =  para$SVM$C, 
	                     kpar = list(sigma = para$SVM$sigma), 
	                     cross = 0, prob.model = TRUE, scaled = TRUE)
	}                   
        pp <-  abs(predict(result$main.model, newdata = X, type = "probabilities"))[,2] 
        result$trn <- pp
        pp <- abs(predict(result$main.model, newdata =X.val[, rhs.vars], type = "probabilities"))[,2] 
        result$val <- pp
        pp <- abs(predict(result$main.model, newdata =X.tst, type = "probabilities"))[,2] 
        result$tst <- pp 
        result$form <- formula.string                 
	return(result)
} 
)
return(res)
}
#' @rdname Train 
#' @export 
Train.Validate.Test <- function(classifier, X.trn, Y.trn, X.val, Y.val, X.tst, Y.tst,  
                           para, opt.para = FALSE, return.model = FALSE, varimp = TRUE, ...){
                           
Train.Models <- lapply(TrainModels(), function(x) x)  

srt <- proc.time()                        
mod <- Train.Models[[classifier]](X.trn, Y.trn, X.val, X.tst, para, opt.para)
ed <- proc.time()
run.time <- as.numeric((ed - srt)[3]) 

collect.garbage()
model <- mod$main.model

Imp <- NULL
if(varimp){
if(classifier %in%c("GBM", "GLM", "RF")){
xx <- varImp(model$finalModel, scale = TRUE)
Imp <- xx[order(xx$Overall, decreasing = TRUE), ,drop = FALSE][ , ,drop = FALSE]
colnames(Imp) <- "importance"
}
}

para =  mod$para
trn.pred <- mod$trn
tst.pred <- mod$tst
val.pred <- mod$val
val.obs <- factor(Y.val[, 1])
tst.obs <-  factor(Y.tst[, 1])

### get optimal threshold using pedicted training probabilities and training performance 
val.perf <- Performance.measures(pred=val.pred, obs=val.obs, prevalence = para$prior)
### get optimal thresholds
thresh <- as.numeric(val.perf$threshold)
para$thresh = thresh 

#### performance on test set using optimal thresholds 
tst.perf <-  Performance.measures(pred=tst.pred, obs=tst.obs, prevalence = para$prior, 
                         threshold=thresh)  
                         
if(return.model)                        
res <- list(model = model, para = mod$para, run.time = run.time, varimp = Imp, pred.prob = tst.pred, true = tst.obs, 
           perf = list(tst.perf=tst.perf))
else 
res <- list(para = mod$para, run.time = run.time, varimp = Imp, perf = list(tst.perf=tst.perf))
return(res)
}























