# nsqipDischargeToRehab
Currently this pacakage is not for public use.

# Example 
```
suppressMessages(require(nsqipDischargeToRehab))
suppressMessages(require(plyr))
suppressMessages(require(caret))
suppressMessages(require(PresenceAbsence))
suppressMessages(require(missForest))
suppressMessages(require(parallel))
suppressMessages(require(doParallel))
cores <- floor(0.8*detectCores()) 
registerDoParallel(cores=cores)
# suppressMessages(require(stargazer))  ## for latex output 

para = list(
  method = "cv", 
  tuneLength = 4,
  number = 3,
  GLM = c(), 
  ELR = list(grid = expand.grid(gamma = c(0.01, 0.25, 10.05), p =  c(100, 400, 800)),  
             ken = "sigmoid", p = 400, gamma = 10.05), 
  GBM = list(n.trees=10, interaction.depth=1, shrinkage=0.01, 
             mstop = 100, n.minobsinnode = 10,prune = "no", 
             grid = expand.grid(interaction.depth = c(1, 3, 5, 10), 
                                n.trees= c(10, 20, 50, 100), shrinkage=c(0.01,0.1), 
                                n.minobsinnode=10)),
  RF = list(ntree = 100, mtry = 5), 
  pLDA = list(lambda = 0.5),
  avNNET = list(nnet.size = 50, nnet.decay= 0.001, nnet.bag = 10, maxit = 200), 
  SVM = list(C = .05, sigma = 0.5, kernel= "rbfdot", 
              grid = expand.grid(C = c(0.01, 1, 10, 100), sigma = c(0.001, 0.05, 0.5))),
  GLMnet = list(alpha = 1, lambda = 1, grid = expand.grid(alpha=c(0,0.5, 1),
                            lambda= c(0.01, 0.05, 0.1, 0.5)))
) 

link <- "nsqip_spine_11_12data.csv"
Data <- preProcessData(file.name =link, impute = FALSE)

XY.dat <- Data$XY.dat
form <- Data$form.list
resp.vars <- Data$resp.vars
rhs.vars <-  Data$rhs.vars 



classifier <- c("GLM", "GBM", "GLMnet", "RF") ## "ELR", "SVM", "avNNET", )
### to run ELR, first intall it using install ELR using 
### devtools::install_github("nguforche/ELR") 
### suppressMessages(require(ELR))
cv = 10 
#ix <- sample(nrow(XY.dat), nrow(XY.dat))
ix <- sample(nrow(XY.dat), 12000)
dat <- XY.dat[ix, ]

Rehab.cv <- RunCV(form=form, data=dat, classifier=classifier, cv = cv, para=para, opt.para=TRUE, 
                  val.propo = 0.10, return.model = FALSE, varimp = TRUE, parallel = TRUE, 
                  mc.cores = cores, seed=12345)



Rehab.para <- get.para(object=Rehab.cv, para=para)
Rehab.perf <-  get.cvResults(object=Rehab.cv)
Rehab.perf
## stargazer(Rehab.perf[, -c(3,4)], summary = FALSE, rownames = FALSE, title = "Performance")
nBoots = 100
Rehab.boot <- RunBootstrap(form=form, data=dat, classifier=classifier,  nBoots = nBoots, para=cv.para, 
                         opt.para=FALSE, val.propo = 0.10, return.model = FALSE, varimp = FALSE, 
                         parallel = TRUE, mc.cores = cores, seed=12345)
boot.perf <-  get.bootResults(object=Rehab.boot)
boot.perf

```
## Install using devtools: 

```sh
devtools::install_github("nguforche/nsqipDischargeToRehab")
```
## Status
Experimental 





