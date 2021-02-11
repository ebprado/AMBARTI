library(devtools)
library(dbarts)
library(AMBARTI)

data = generate_data_AMMI(10, 10, 1,1,1,c(8, 10, 12))
db = data.frame(y = data$y, data$x)
x = data$x
y = data$y
mod1 = lm(y ~ g + e, data=db)
res = mod1$residuals

xx <- model.matrix(~ -1 + g + e, data=x,
                      contrasts.arg=list(g=contrasts(x$g, contrasts=F),
                                         e=contrasts(x$e, contrasts=F)))

# write.csv(xx, 'design_matrix.csv')
### aux g
x_aux_g <- model.matrix(~ (-1 + .)^2, data=data.frame(xx[, grepl('g', colnames(xx))]))[,-c(1:10)]
x_g = matrix(NA, ncol=ncol(x_aux_g), nrow=nrow(x_aux_g))
colnames(x_g) = colnames(x_aux_g)
for (k in 1:ncol(x_g)){
  name_col_g = unlist(strsplit(colnames(x_g)[k],':'))
  x_g[,k] = apply(xx[,name_col_g],1,sum)*3
}


### aux e
x_aux_e <- model.matrix(~ (-1 + .)^2, data=data.frame(xx[, grepl('e', colnames(xx))]))[,-c(1:10)]
x_e = matrix(NA, ncol=ncol(x_aux_e), nrow=nrow(x_aux_e))
colnames(x_e) = colnames(x_aux_e)
for (k in 1:ncol(x_e)){
  name_col_e = unlist(strsplit(colnames(x_e)[k],':'))
  x_e[,k] = apply(xx[,name_col_e],1,sum)*(-2)
}


## aux x_g_e
x_aux_g_e = cbind(x_g, x_e)
x_g_e = matrix(NA, ncol=ncol(x_g) * ncol(x_e), nrow = nrow(x_g))
colnames(x_g_e) = paste('A', seq(1:ncol(x_g_e)))
k=0
for (i in 1:ncol(x_g)){
  for (j in 1:ncol(x_e)){
    k = k+1
    names_col = c(colnames(x_g)[i], colnames(x_e)[j])
    x_g_e[,k] =  apply(x_aux_g_e[,names_col], 1, sum)
    colnames(x_g_e)[k] = paste(names_col, collapse = ':')
  }
}

dim(x_g_e)

test = as.data.frame(x_g_e[, sample(1:ncol(x_g_e), size = 1)])

bartfit = dbarts::bart2(x_g_e,res, n.chains = 4, n.burn = 10000, n.samples = 1000, sigest = 1, base = 0.9, power = 0.3)

cor(y, mod1$fitted.values + bartfit$yhat.train.mean)


ambart = AMBARTI::ambarti(x, y, ntrees = 50, nburn = 100, npost = 100)
ambart_yhat = apply(ambart$y_hat,2,mean)
apply(ambart$beta_hat,2,mean)
mod1$coefficients
cor(y, ambart_yhat)
cor(y, mod1$fitted.values)
