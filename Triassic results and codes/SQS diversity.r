# created by Xioakang Liu on 06/05/2025 => xkliu@cug.edu.cn
#################################################
####################ostracods####################
#################################################
library(divDyn)
g.dat <- read.csv("D:/PyRate-master/DeepMorphT/ostracods.csv")
g.dat$stg <- as.integer(g.dat$stg)
g.dat <- na.omit(g.dat)
#data(stages)
stages <- read.csv("stages.csv")
# raw diversity
g.div <- divDyn(g.dat, tax = "genus", bin = "stg")
# SQS sampling
p.sqs.div <- subsample(x = g.dat,iter = 1000,q = 0.7,tax = "genus",bin = "stg",type = "sqs",output = "dist")#,keep=c(49,58)
subArit <- subsample(g.dat, bin="stg", tax="genus", iter=1000, q=0.7, output="arit",type = "sqs")

# results
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Ost Range-through diversity (genera)", ylim = c(0, 150))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(p.sqs.div$divRT[51:55,]), col = "orange", lwd = 2)
lines(stages$mid[51:55], subArit$divRT[51:55], col = "blue", lwd = 2)
shades(stages$mid[1:58], p.sqs.div$divRT, col="orange", res=c(0.05,0.25,0.75,0.95))
legend("topright", legend = c("Raw", "SQS (q=0.7)","subArit"),
       col = c("black", "orange", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Ost gap filler Rates", ylim = c(0, max(c(p.sqs.div$oriGF, p.sqs.div$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], p.sqs.div$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], p.sqs.div$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("Origination Rate", "Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")
#write.csv(p.sqs.div$extGF[51:55,], "p.sqs.div$extGF[4958,].csv")#, row.names = FALSE)

# cr diversity
cr40 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=40, duplicates=TRUE, iter=1000)
cr80 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=80, duplicates=TRUE, iter=1000)
cr80arit <- subsample(g.dat, tax="genus", bin="stg", output = "arit",q=80, duplicates=TRUE, iter=1000)
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Ost Range-through diversity (genera)", ylim = c(0, 150))
shades(stages$mid[1:58], cr80$divRT, col="blue", res=c(0.05,0.25,0.75,0.95))
shades(stages$mid[1:58], cr40$divRT, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(cr40$divRT[51:55,]), col = "red", lwd = 2)
lines(stages$mid[51:55], rowMeans(cr80$divRT[51:55,]), col = "blue", lwd = 2)
legend("topright", legend = c("Raw", "CR40", "CR80"),
       col = c("black", "red", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Ost  gap filler Rates", ylim = c(0, max(c(cr80$oriGF, cr80$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], cr80$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], cr80$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("CR80 Origination Rate", "CR80 Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")

write.csv(cr80$extGF[51:55,], "cr80$extGF[4958,].csv")#, row.names = FALSE)
write.csv(rarefDD$extGF[51:55,], "cr80$extGF[4958,].csv")#, row.names = FALSE)

#################################################
###################brachiopods###################
#################################################
library(divDyn)
g.dat <- read.csv("D:/PyRate-master/DeepMorphT/brachiopods.csv")
g.dat$stg <- as.integer(g.dat$stg)
g.dat <- na.omit(g.dat)
# data(stages)
stages <- read.csv("stages.csv")
# raw diversity
g.div <- divDyn(g.dat, tax = "genus", bin = "stg")
# SQS sampling
p.sqs.div <- subsample(x = g.dat,iter = 1000,q = 0.7,tax = "genus",bin = "stg",type = "sqs",output = "dist")#,keep=c(49,58)
subArit <- subsample(g.dat, bin="stg", tax="genus", iter=1000, q=0.7, output="arit",type = "sqs")
# results
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Bra Range-through diversity (genera)", ylim = c(0, 350))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(p.sqs.div$divRT[51:55,]), col = "orange", lwd = 2)
lines(stages$mid[51:55], subArit$divRT[51:55], col = "blue", lwd = 2)
shades(stages$mid[1:58], p.sqs.div$divRT, col="orange", res=c(0.05,0.25,0.75,0.95))
legend("topright", legend = c("Raw", "SQS (q=0.7)","subArit"),
       col = c("black", "orange", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Bra gap filler Rates", ylim = c(0, max(c(p.sqs.div$oriGF, p.sqs.div$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], p.sqs.div$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], p.sqs.div$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("Origination Rate", "Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")
#write.csv(p.sqs.div$extGF[51:55,], "p.sqs.div$extGF[4958,].csv")#, row.names = FALSE)

# cr diversity
cr40 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=40, duplicates=TRUE, iter=1000)
cr80 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=80, duplicates=TRUE, iter=1000)
cr80arit <- subsample(g.dat, tax="genus", bin="stg", output = "arit",q=80, duplicates=TRUE, iter=1000)
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Bra Range-through diversity (genera)", ylim = c(0, 350))
shades(stages$mid[1:58], cr80$divRT, col="blue", res=c(0.05,0.25,0.75,0.95))
shades(stages$mid[1:58], cr40$divRT, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(cr40$divRT[51:55,]), col = "red", lwd = 2)
lines(stages$mid[51:55], rowMeans(cr80$divRT[51:55,]), col = "blue", lwd = 2)
legend("topright", legend = c("Raw", "CR40", "CR80"),
       col = c("black", "red", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Bra  gap filler Rates", ylim = c(0, max(c(cr80$oriGF, cr80$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], cr80$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], cr80$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("CR80 Origination Rate", "CR80 Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")

# write.csv(cr80$extGF[51:55,], "cr80$extGF[4958,].csv")#, row.names = FALSE)
# write.csv(rarefDD$extGF[51:55,], "cr80$extGF[4958,].csv")#, row.names = FALSE)
# 

#################################################
####################ammonoids####################
#################################################
library(divDyn)
g.dat <- read.csv("D:/PyRate-master/DeepMorphT/ammonoids.csv")
g.dat$stg <- as.integer(g.dat$stg)
g.dat <- na.omit(g.dat)
# data(stages)
stages <- read.csv("stages.csv")
# raw diversity
g.div <- divDyn(g.dat, tax = "genus", bin = "stg")
# SQS sampling
#p.sqs.div <- subsample(x = g.dat,iter = 1000,q = 0.7,tax = "genus",bin = "stg",type = "sqs",output = "dist")#,keep=c(49,58)
subArit <- subsample(g.dat, bin="stg", tax="genus", iter=1000, q=0.7, output="arit",type = "sqs")
# results
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Amm Range-through diversity (genera)", ylim = c(0, 300))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(p.sqs.div$divRT[51:55,]), col = "orange", lwd = 2)
lines(stages$mid[51:55], subArit$divRT[51:55], col = "blue", lwd = 2)
shades(stages$mid[1:58], p.sqs.div$divRT, col="orange", res=c(0.05,0.25,0.75,0.95))
legend("topright", legend = c("Raw", "SQS (q=0.7)","subArit"),
       col = c("black", "orange", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Amm gap filler Rates", ylim = c(0, max(c(p.sqs.div$oriGF, p.sqs.div$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], p.sqs.div$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], p.sqs.div$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], subArit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("Origination Rate", "Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")
#write.csv(p.sqs.div$extGF[51:55,], "p.sqs.div$extGF[4958,].csv")#, row.names = FALSE)

# cr diversity
cr40 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=40, duplicates=TRUE, iter=1000)
cr80 <- subsample(g.dat, tax="genus", bin="stg", output = "dist",q=80, duplicates=TRUE, iter=1000)
cr80arit <- subsample(g.dat, tax="genus", bin="stg", output = "arit",q=80, duplicates=TRUE, iter=1000)
tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Amm Range-through diversity (genera)", ylim = c(0, 300))
shades(stages$mid[1:58], cr80$divRT, col="blue", res=c(0.05,0.25,0.75,0.95))
shades(stages$mid[1:58], cr40$divRT, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], g.div$divRT[51:55], lwd = 2)
lines(stages$mid[51:55], rowMeans(cr40$divRT[51:55,]), col = "red", lwd = 2)
lines(stages$mid[51:55], rowMeans(cr80$divRT[51:55,]), col = "blue", lwd = 2)
legend("topright", legend = c("Raw", "CR40", "CR80"),
       col = c("black", "red", "blue"), lwd = c(2, 2), bg = "white")

tsplot(stages, shading = "stage", boxes = "sys", xlim = c(254.14, 237),
       ylab = "Amm  gap filler Rates", ylim = c(0, max(c(cr80$oriGF, cr80$extGF), na.rm = TRUE)))
# origination rates
shades(stages$mid[1:58], cr80$oriGF, col="blue", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$oriGF[51:55], col = "blue", lwd = 2)
# extinction rates
shades(stages$mid[1:58], cr80$extGF, col="red", res=c(0.05,0.25,0.75,0.95))
lines(stages$mid[51:55], cr80arit$extGF[51:55], col = "red", lwd = 2)
legend("topright", legend = c("CR80 Origination Rate", "CR80 Extinction Rate"),
       col = c("blue", "red"), lwd = 2, bg = "white")