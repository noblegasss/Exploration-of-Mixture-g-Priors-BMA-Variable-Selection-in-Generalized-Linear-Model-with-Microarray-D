## Load package
library(genefilter)
library(annotate)
library(MASS)
library(BAS)
library(caret)

# Prepare data
## Load data
prostate = read.csv("Final Project/Prostate_GSE6919_U95C.csv", row.names = 1, header = F)
colnames(prostate) = prostate[1,]
prostate = prostate[-1,]

## Filter probes and annotation
exp = na.omit(t(prostate))[-1,]
arrayIQR<-apply(exp,1,IQR)
probe<-rownames(exp)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,'hgu95c.db')
exp2<-exp[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu95c.db")
rownames(exp2)<-geneSymbol
### Scale 
exp22 = apply(exp2, 1, as.numeric)
exp22 = scale(exp22)
rownames(exp22) = colnames(exp2)
prostate.new = data.frame(type = prostate$type, exp22)
prostate.new$type = ifelse(prostate.new$type == 'primary_prostate_tumor', 1, 0)

## Select top 10% variation genes
prostate.new2 = varFilter(as.matrix(t(prostate.new[,-1])),var.cutoff = 0.9)
prostate.new2 = data.frame(type = prostate.new$type, t(prostate.new2))

# Analysis
## Run for 1 sample 
sam = sample.int(115, 92)
train = prostate.new2[sam,]
test = prostate.new2[-sam,]
## Fit BAS using CCH prior
bas1 = bas.glm(type~., data = train, family  = binomial, 
               betaprior = beta.p[[1]],
               method = "MCMC+BAS", n.models = 2^7, laplace = T)
pred = predict(bas1, test, type = "BPM") ### Best predictor model
pred0 = ifelse(1/(1+exp(-pred$fit))>=0.5,1,0)
mean(pred0 != test$type) ### classification error 
mean((1/(1+exp(-pred$fit)) - test$type)^2) ### Brier score 
image(bas1, rotate = F)
s = summary(bas1)
s01 = na.omit(s[s[,1] > 0.01,1])


## CV for different gpriors and model prior
logit = function(x) 1/(1+exp(-x))
cv.bas = function(data, beta.p, model.p, k = 10, keep.prob = 0.01){
  fold = createFolds(data[,1], k = k)
  model = lapply(fold, function(i) bas.glm(type~., data = data[-i,], family = binomial,
                                           betaprior = beta.p, modelprior = model.p,
                                           method = "MCMC+BAS", n.models = 2^7, laplace = T))
  b.select = lapply(1:k, function(j){
    s0 = summary(model[[j]])
    names(na.omit(s0[s0[,1] > keep.prob,1]))})
  pred = sapply(1:k, function(d){
    pred0 = predict(model[[d]], data[fold[[d]],], type = "BPM")$fit
    label = ifelse(logit(pred0) >= 0.5, 1, 0)
    c(mean((logit(pred0) - data[fold[[d]],1])^2, na.rm = T),
    mean(label != data[fold[[d]],1], na.rm = T))
  })
  cv.error = rowMeans(pred)
  names(cv.error) = c("brier score", "test error")
  select = table(unlist(b.select))
  return(list(model = model, b.select = select, cv.error = cv.error))
}

## Beta binomial(1,1) prior
model.p = beta.binomial(1,1)
n = as.numeric(nrow(prostate.new2))
p = as.numeric(ncol(prostate.new2))
beta.p = list(CCH(0.5, n, s = 0), 
              CCH(0.5, p, s = 0),
              CCH(1, n, s = 0),
              CCH(1, p, s = 0),
              robust(),
              hyper.g(),
              hyper.g.n(),
              TG(),
              testBF.prior(p),
              intrinsic(),
              IC.prior(n),
              IC.prior(2))

model.all = NULL
b.select.all = list()
cv = matrix(nrow = 12, ncol = 2)
for(i in 1:12){
  l = cv.bas(data=prostate.new2, beta.p = beta.p[[i]], model.p = model.p)
  model.all = append(model.all, l$model) 
  b.select.all[[i]] = l$b.select
  cv[i,] = l$cv.error
}
  
beta.all = unique(names(unlist(b.select.all)))
b.all0 = matrix(0, nrow = length(beta.all), ncol = 12)
rownames(b.all0) = beta.all

for(i in 1:12){
  b.all = c(b.all0[,i], b.select.all[[i]])
  b.all01 = tapply(b.all, names(b.all), sum)
  b.all0[,i] = b.all01[match(beta.all,names(b.all01) )]
}

colnames(b.all0) = c("CCH0.5n", "CCH0.5p", "CCH1n", "CCH1p", "robust", "hyper.g",
                     "hyper.g.n", "TG","TBF", "intrinsic", "BIC", "AIC")

b.all0 = b.all0[rownames(b.all0)!= "Intercept", ]
dat <- data.frame(b.all0) %>%
  rownames_to_column('Seleted.genes') %>%
  gather(gpriors, value, -Seleted.genes) 

ggplot(dat, aes(x=gpriors,y=Seleted.genes,fill=value)) + geom_tile()+ 
  scale_fill_gradient(low = "white", high = "navy") +
  theme(axis.text.y = element_text(size=5))


## uniform prior
model.p2 = uniform()
n = as.numeric(nrow(prostate.new2))
p = as.numeric(ncol(prostate.new2))

model.all2 = NULL
b.select.all2 = list()
cv2 = matrix(nrow = 12, ncol = 2)
for(i in 1:12){
  l = cv.bas(data=prostate.new2, beta.p = beta.p[[i]], model.p = model.p2)
  model.all2 = append(model.all, l$model) 
  b.select.all2[[i]] = l$b.select
  cv2[i,] = l$cv.error
}

beta.all2 = unique(names(unlist(b.select.all2)))
b.all02 = matrix(0, nrow = length(beta.all2), ncol = 12)
rownames(b.all02) = beta.all2

for(i in 1:12){
  b.all = c(b.all02[,i], b.select.all2[[i]])
  b.all01 = tapply(b.all, names(b.all), sum)
  b.all02[,i] = b.all01[match(beta.all2,names(b.all01) )]
}

colnames(b.all02) = c("CCH0.5n", "CCH0.5p", "CCH1n", "CCH1p", "robust", "hyper.g",
                     "hyper.g.n", "TG","TBF", "intrinsic", "BIC", "AIC")

b.all02 = b.all02[rownames(b.all02)!= "Intercept", ]
dat2 <- data.frame(b.all02) %>%
  rownames_to_column('Seleted.genes') %>%
  gather(gpriors, value, -Seleted.genes) 

ggplot(dat2, aes(x=gpriors,y=Seleted.genes,fill=value)) + geom_tile()+ 
  scale_fill_gradient(low = "white", high = "navy") + 
  theme(axis.text.y = element_blank()) 


## Summary

unique.b = sapply(1:12, function(b) length(b.select.all[[b]]))
unique.b2 = sapply(1:12, function(b) length(b.select.all2[[b]]))

dat.bebi = rbind(unique.b, t(cv))
colnames(dat.bebi) = c("CCH0.5n", "CCH0.5p", "CCH1n", "CCH1p", "robust", "hyper.g",
                       "hyper.g.n", "TG","TBF", "intrinsic", "BIC", "AIC")
rownames(dat.bebi) = c("TotalNumberOfB", "cv.error", "cv.brier")
xtable(dat.bebi)

dat.bebi2 = rbind(unique.b2, t(cv2))
colnames(dat.bebi2) = c("CCH0.5n", "CCH0.5p", "CCH1n", "CCH1p", "robust", "hyper.g",
                       "hyper.g.n", "TG","TBF", "intrinsic", "BIC", "AIC")
rownames(dat.bebi2) = c("TotalNumberOfB", "cv.brier", "cv.error")
xtable(dat.bebi2)

par(mfrow = c(1,2))
image(model.all[[9]],rotate = F)
image(model.all[[7]],rotate = F)


## Compare with elastic net and lasso

model.enet = function(X, y, k = 10){
  fold = createFolds(data[,1], k = k)
  a <- seq(0.1, 0.9, 0.05)
  model = lapply(fold, function(j){
    search <- foreach(i = a, .combine = rbind) %dopar% {
    cv <- cv.glmnet(X[-j,], y[-j], family = "binomial", nfold = 10,paralle = TRUE, alpha = i)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
    }
  cv <- search[search$cvm == min(search$cvm), ]
  glmnet(X[-j,], y[-j], family = "binomial", lambda = cv$lambda.min, alpha = cv$alpha)
  })
  b.select = lapply(1:k, function(i){
    w = data.frame(as.matrix(coef(model[[i]])))
    b = w %>% filter(s0!=0)
    rownames(b)
    })
  pred = sapply(1:k, function(d){
    pred0 = predict(model[[d]], X[fold[[d]],], type = "response")
    label = ifelse(logit(pred0) >= 0.5, 1, 0)
    c(mean((logit(pred0) - y[fold[[d]]])^2, na.rm = T),
      mean(label != y[fold[[d]]], na.rm = T))
  })
  cv.error = rowMeans(pred)
  names(cv.error) = c("brier score", "test error")
  select = table(unlist(b.select))
  return(list(model = model, b.select = b.select, unique.select = select, cv.error = cv.error))
}

enet = model.enet(X = as.matrix(prostate.new2[,-1]), y = prostate.new2[,1])

model.lasso = function(X, y, k = 10){
  fold = createFolds(data[,1], k = k)
  model = lapply(fold, function(j){
    cv <- cv.glmnet(X[-j,], y[-j], family = "binomial", nfold = 10,alpha = 1)
    glmnet(X[-j,], y[-j], family = "binomial", lambda = cv$lambda.min, alpha = 1)
  })
  b.select = lapply(1:k, function(i){
    w = data.frame(as.matrix(coef(model[[i]])))
    b = w %>% filter(s0!=0)
    rownames(b)
  })
  pred = sapply(1:k, function(d){
    pred0 = predict(model[[d]], X[fold[[d]],], type = "response")
    label = ifelse(logit(pred0) >= 0.5, 1, 0)
    c(mean((logit(pred0) - y[fold[[d]]])^2, na.rm = T),
      mean(label != y[fold[[d]]], na.rm = T))
  })
  cv.error = rowMeans(pred)
  names(cv.error) = c("brier score", "test error")
  select = table(unlist(b.select))
  return(list(model = model, b.select = b.select, unique.select = select, cv.error = cv.error))
}

lasso = model.lasso(X = as.matrix(prostate.new2[,-1]), y = prostate.new2[,1])

dat.all = data.frame(dat.bebi[,c(1,7,9,10)], lasso = c(length(lasso$unique.select), lasso$cv.error), enet = c(length(lasso$unique.select), lasso$cv.error))
xtable(dat.all)





