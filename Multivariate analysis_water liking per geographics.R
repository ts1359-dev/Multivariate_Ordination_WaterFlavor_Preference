rm(list = ls())
# 1.Ordination analysis of mean water liking (i.e., flavor preference).
# Variables are continuous. Principal component analysis (PCA) is the most appropriate ordination analysis for this data.
library(readxl)
DATA <- as.data.frame(read_excel("C:/Users/vt135/Downloads/Data analysis for env engr/Assignment5/Data_for_Assignment5.xlsx"))
library(robustHD)
ZDATA = standardize(DATA,centerFun = mean,scaleFun = sd)
n = nrow(ZDATA)
Z = as.matrix(ZDATA)
A = (1/(n-1))*crossprod(Z)
e <- eigen(A)
lam_ordered = e$values
lam_ordered
#eigen values: num [1:5] 3.573808  1.306574  1.196178e-01  7.071899e-16  2.775558e-17
library(pracma)
k=ncol(ZDATA)

X = matrix(NA,k,k)
for (i in 1:k){
  ALI = A-(lam_ordered[i]*diag(k))
  X[,i] = null(ALI)
}
X = as.matrix(e$vectors)
X
#eigen vectors are shown below:


[,1]       [,2]         [,3]        [,4]          [,5]
[1,] 0.4778847 -0.2897941  0.787115006  0.26095137  0.000000e+00
[2,] 0.4777971 -0.3293365 -0.595513296  0.55552796 -3.944652e-15
[3,] 0.4902821 -0.3249549 -0.157622382 -0.78329288 -1.250409e-01
[4,] 0.3593924  0.6419235 -0.004093671  0.06706119 -6.739861e-01
[5,] 0.4168894  0.5384193 -0.030859492 -0.07244397  7.280848e-01

PC_SCORES = as.matrix(ZDATA)%*%X
PC_SCORES
[,1]        [,2]         [,3]          [,4]          [,5]
[1,] -1.3370517  0.36965421 -0.323559039 -2.567391e-16 -4.440892e-16
[2,] -0.9043436 -0.12463099  0.085052603 -2.081668e-17  2.220446e-16
[3,] -0.1093546 -0.56948898  0.121247636 -3.642919e-16  1.110223e-16
[4,] -0.9140244 -0.06920544 -0.059324872  7.979728e-17  0.000000e+00
[5,]  0.2191966 -1.77279297 -0.009267016 -1.249001e-16  0.000000e+00
[6,]  1.2821967 -1.80758716  0.015183440 -9.020562e-17 -1.110223e-16
[7,] -0.2385174 -0.39659289 -0.077999961  4.857226e-17 -1.942890e-16
[8,]  1.4921183 -1.87625156 -0.308208738 -2.220446e-16 -4.996004e-16
[9,] -2.4743187 -0.26495953 -0.037393522 -6.800116e-16  1.110223e-16
[10,]  1.2263342 -0.39777456  0.568153628 -5.793976e-16  8.049117e-16
[11,]  1.3989073  0.62663500  0.669876743 -6.800116e-16  9.992007e-16
[12,]  1.6505441  1.34350121 -0.584259606  2.359224e-16 -8.881784e-16
[13,]  2.0219692 -0.95054442 -0.289079594 -2.220446e-16 -7.771561e-16
[14,]  2.8238676  1.10086473  0.491316093 -8.326673e-16  4.440892e-16
[15,]  3.0256020  1.18560153 -0.453271971 -4.440892e-16 -8.881784e-16
[16,]  0.6024343  2.59608040  0.001265626 -3.053113e-16  0.000000e+00
[17,] -2.8894795  0.48009561  0.488973700 -5.689893e-16  1.110223e-15
[18,] -3.0148075  1.05483052 -0.247300876 -7.632783e-17 -1.665335e-16
[19,] -1.0103147 -0.47202873  0.192150224 -7.702172e-16  4.996004e-16
[20,] -2.8509580 -0.05540600 -0.243554500 -1.665335e-16 -2.220446e-16

perVE = lam_ordered/sum(lam_ordered)
perVE
# PC1 explains 71.5% of variance.
# PC2 explains 26.1% of variance.

PC_SCORES = as.matrix(ZDATA)%*%X
PC_SCORES
Xnames = c('MeanLikingSp', 'MeanLikingFr', 'MeanLikingIt', 
           'MeanLikingPor', 'MeanLikingGr')

#ylabs = vector labels
biplot(PC_SCORES[,1:2],X[,1:2],ylabs = Xnames,
       xlab = "PC1", ylab = "PC2")
# In the plot, on the x-axis, positive values mean that people like their type of water whereas negative values mean that people dislike 
#the type of water serviced to them.
# In the plot, on the y-axis, positive values mean that people like the water types in Portugal and Greece and do not like the water types 
#in Spanish, France, and Italy as much. Conversely, negative values imply that people like the water types in Spanish, France, and Italy 
#and do not like the water types in Portugal and Greece as much.

#2.Perform constrained ordination. Build an MLR on each mean liking variable using all of the variables included in the chemistry data that do not exhibit 
#multicollinearity as independent variables. If any chemistry variables violate the VIF rule remove the largest value and re-evaluate 
#until all values are compliant. From your MLR analysis you should get a set of z + 1 (for the intercept) beta coefficients for each mean 
#liking variable. Use these beta coefficients and your chemistry data to predict new values for each mean liking variable. 
#Report the first 5 x 5 segment of the matrix of the predicted liking values

library(FactoMineR)
library(factoextra)

res.pca = PCA(ZDATA, graph = FALSE,ncp = 13)
fviz_pca_biplot(res.pca,repel = TRUE,label = "var",
                col.var = "black")
#TASK2
CONSTRAINTS = as.data.frame(read_excel("C:/Users/vt135/Downloads/Data analysis for env engr/Assignment5/Data_for_Assignment5.xlsx",sheet = "Chemistry"))

CDATA = standardize(CONSTRAINTS,centerFun = mean,
                    scaleFun = sd)
n = nrow(CDATA) 
COR_C = (1/(n-1))*(t(CDATA)%*%as.matrix(CDATA))
INV = Matpow(COR_C,-1) 
VIF = diag(INV)
VIF
#remove TDS and repeat
CDATA2 = as.matrix(CDATA[,2:6])
COR_C = (1/(n-1))*(t(CDATA2)%*%CDATA2)
INV = Matpow(COR_C,-1)  
VIF = diag(INV)
VIF                   
u = ncol(CDATA2)+1
B = matrix(NA,u,ncol(ZDATA))

MASTER = as.data.frame(cbind(ZDATA,CDATA2))
for (i in 1:ncol(ZDATA)){
  FIT = lm(MASTER[,i] ~ MASTER$Cl+MASTER$SO42+MASTER$HCO3+
             MASTER$Na+ + MASTER$pH)
  B[,i] = FIT$coefficients
}
PREDICTED = matrix(NA,n,ncol(ZDATA))
for (i in 1:ncol(ZDATA)){
  PREDICTED[,i] = B[1,i]+ B[2,i]*CDATA2[,1] + 
    B[3,i]*CDATA2[,2] + B[4,i]*CDATA2[,3]+ B[5,i]*CDATA2[,4]+B[6,i]*CDATA2[,5]
}
print(PREDICTED[1:5, 1:5])
PREDICTED = as.data.frame(PREDICTED)
CDATA_PRED = standardize(PREDICTED)
cnames = c('MeanLikingSp','MeanLikingFr','MeanLikingIt','MeanLikingPor', 'MeanLikingGr')
colnames(CDATA_PRED)<- cnames
PCfit = PCA(CDATA_PRED,graph = FALSE,ncp = 5)
EIGVEC = PCfit$svd$V
EIGVEC
SCORES = as.matrix(CDATA_PRED)%*%as.matrix(EIGVEC) 
SCORES

v = nrow(CDATA_PRED)
w = as.matrix(CDATA_PRED)
A = (1/(v-1))*crossprod(w)
e <- eigen(A)
lam_ord = e$values
percV = lam_ord/sum(lam_ord)
percV
# The 1st mode explains 74.4% of the predicted variance.
# The 2nd mode explains 24.7% of the predicted variance.

L = PCfit$eig[,1]
Cmat = matrix(NA,ncol(CDATA2),2)

for (i in 1:(u-1)){
  c1 = sqrt(L[1]/sum(L)) #VE Correction for PC1
  c2 = sqrt(L[2]/sum(L)) #VE Correction for PC2
  Cmat[i,1:2] = c(cor(CDATA2[,i],SCORES[,1])*c1, 
                  cor(CDATA2[,i],SCORES[,2])*c2)
}
cnames2 = c('Cl','SO42','HCO3','Na+', 'pH')
Cmat = as.data.frame(Cmat)
rownames(Cmat)<- cnames2
p2=fviz_pca_biplot(PCfit,repel = TRUE)
fviz_add(p2,as.data.frame(Cmat*4.5),geom = "arrow",
         repel = TRUE)
#Interpretation of PCA_constrained ordination biplot:
#1.A significant number of people (in the lower left quadrant) do not like the type of water in any of the countries. The type of water 
#they prefer has low pH, low concentration of SO42- and moderately low concentration of Cl, Na+, and HCO3.

#2.Some people (in the upper left quadrant) do not like the type of water in any of the countries. The type of water they prefer has low pH,
#low concentration in HCO3 and SO42- and high concentration of Cl, Na+.

#3.A significant number of people (in the lower right quadrant) like the type of water in France, Italy, and Spain, which has 
#low concentration of Cl, Na+, and SO42-, is slightly high in concentration of HCO3 -, and has high pH. This group of people does not like
#the type of water in Portugal or Greece.

#4.Some people (in the upper right quadrant) like the type of water in Portugal and Greece, which has low concentration of Cl, Na+, and 
#HCO3, high concentration of SO42-, and has low pH. This group of people does not like the type of water in France, Italy, or Spain.

