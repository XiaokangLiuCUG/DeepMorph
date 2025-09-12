# created by Xioakang Liu on 10/07/2024 => xkliu@cug.edu.cn
library(geomorph)
library(dispRity)
library(factoextra)
#taxa <- as.list(read.table("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/ammonoid/taxa.txt",sep='\t'))

df = read.csv("C:/Users/xiaokang/Jupyter_code/PTDisparity_recovery/landmarks/ammonoid/ammonoid_a_landmarks_alldx.csv",header=F,row.names = 1)
ammonoid.links <- as.matrix(read.table("C:/Users/xiaokang/Jupyter_code/PTimages/landmarks/brachiopod/brachiopod_links.txt"))
curves<-as.matrix(read.csv("C:/Users/xiaokang/Jupyter_code/PTimages/landmarks/brachiopod/brachiopod_curveslide.csv",header = T))
num_sp = length(row.names(df))
AA <- array(data = unlist(t(df[,1:128])),dim = c(2,64,num_sp),dimnames=NULL)

# df = read.csv("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/Ammonoidea_landmarks2.csv",header=T,row.names = 1)
# ammonoid.links <- as.matrix(read.table("links.txt"))
# curves<-as.matrix(read.csv("amm_curveslide.csv",header = T))
# AA <- array(data = unlist(t(df[,1:84])),dim = c(2,42,113),dimnames=NULL)


BB <- aperm(AA,c(2,1,3))

#BB <- array(BB,dimnames=list(c(1:42),c('x','y'),list(df$V1)))
X.gpa <- gpagen(BB, curves=curves,PrinAxes = FALSE,Proj=T)#PrinAxes = FALSE,Proj=T, approxBE=T)
#X.gpa <- gpagen(AA,PrinAxes = FALSE,ProcD = TRUE)
summary(X.gpa)
plot(X.gpa)

##save the procrustes coordinates to csv
# coo <- X.gpa$coords
# COO2 <- aperm(coo,c(2,1,3))
# COO3 <- array(unlist(COO2),dim = c(113,84))
# #COO3 <- aperm(COO3,c(2,1))
# write.csv(X.gpa$coords,"C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/coo2.csv")

average.shape <- mshape(X.gpa$coords)#average shape points
# globalIntegration(X.gpa$coords)

#lab=c(1:84)

PCA_data<-gm.prcomp(X.gpa$coords)
#write.csv(PCA_data$x[,1:2],"C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/ammonoid/pca_coo_amm.csv")
#PCA_data$x#A
#plot(X.gpa)
plot(PCA_data)
#text(X,Y,labels = lab,pos=4)
#fviz_eig(PCA_data, addlabels = T)
summary(PCA_data)
plotRefToTarget(PCA_data$shapes$shapes.comp1$max,average.shape)#, links = ostracod.links,label=FALSE
plotAllSpecimens(X.gpa$coords)#, links = ammonoid.links

#PCA <- plot.gm.prcomp(PCA_data$coords, axis1 = 1, axis2 = 2, warpgrids = TRUE, label = TRUE)# plot.gm.prcomp plotTangentSpace
#PCA_data$pc.summary

pvar <- (PCA_data$sdev^2)/(sum(PCA_data$sdev^2))*100
names(pvar) <- seq(1:length(pvar))
par(mgp=c(2,1,0))
bar_centers <- barplot(pvar[1:10], col=rgb(0.2,0.4,0.6,0.6), ylim = c(0, max(pvar[1:10]) + 5), xlab= "Principal Components", ylab = "% Variance") # save this plot
text(bar_centers, pvar[1:10] + 2, labels = round(pvar[1:10],2),cex = 0.8)
title("Proportion of variance", line = 1)#, row=2
pvar[1:10]

Group <- read.csv("C:/Users/xiaokang/Jupyter_code/PTDisparity_recovery/landmarks/ammonoid/ammonoid_interval_all-dx2.csv",row.names = 1)
families <- Group[,'Family']

changhsingian <- row.names(Group)[Group[,'Changhsingian'] == "1"]
PTTB <- row.names(Group)[Group[,'PTTB'] == "1"]
LGriesbachian <- row.names(Group)[Group[,'LGriesbachian'] == "1"]
Dienerian <- row.names(Group)[Group[,'Dienerian'] == "1"]

Induan <- row.names(Group)[Group[,'Induan'] == "1"]
Smithian <- row.names(Group)[Group[,'Smithian'] == "1"]
Spathian <- row.names(Group)[Group[,'Spathian'] == "1"]
EAnisian <- row.names(Group)[Group[,'EAnisian'] == "1"]
LAnisian <- row.names(Group)[Group[,'LAnisian'] == "1"]
ELadinian <- row.names(Group)[Group[,'ELadinian'] == "1"]
LLadinian <- row.names(Group)[Group[,'LLadinian'] == "1"]

Olenekian <- row.names(Group)[Group[,'Olenekian'] == "1"]
Anisian <- row.names(Group)[Group[,'Anisian'] == "1"]
Ladinian <- row.names(Group)[Group[,'Ladinian'] == "1"]

plot(PCA_data$x[,1], PCA_data$x[,2],  xlab = "PC1", ylab = "PC2", pch = 20, cex=1)
text(x = PCA_data$x[,1], y = PCA_data$x[,2],labels = c(1:num_sp),cex = 0.8,col = 'purple',srt = 45)
abline(v=c(0), lwd=0.7, col="gray40", lty=2) # add central gridlines 
abline(h=c(0), lwd=0.7, col="gray40",lty=2) # add central gridlines
#title("all species")

cat(max(PCA_data$x[,1]),min(PCA_data$x[,1]))
cat(max(PCA_data$x[,2]),min(PCA_data$x[,2]))


Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor,interval, pch,D){
  par(new=D)
  plot(xcoord, ycoord,  xlab = "PC1", ylab = "PC2",pch = pch, cex=2,axes = !D, xlim = c(min(PCA_data$x[,1]),max(PCA_data$x[,1])), ylim = c(min(PCA_data$x[,2]),max(PCA_data$x[,2])))
  points(mean(xcoord),mean(ycoord), pch = pch+15, cex = 2, xlim = c(-0.3,0.5), ylim = c(min(PCA_data$x[,2]),max(PCA_data$x[,2])))
  cat(mean(xcoord),mean(ycoord))
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.4) , border = 0)
  #title(interval)
  
}
#Plot_ConvexHull(xcoord= PCA_data$x[,1], ycoord= PCA_data$x[,2], lcolor="#FF6666", bgcolor="#FF6666","ALL",pch=0,D=FALSE)

# plot(PCA_data$x[changhsingian,][,1], PCA_data$x[changhsingian,][,2],  xlab = "PC1", ylab = "PC2", xlim = c(-0.42,0.3), ylim = c(-0.35,0.2))#,pch = 20, cex=0.5
# #text(x = PCA_data$x[,1], y = PCA_data$x[,2],labels = seq(1:num_sp),cex = 1,col = 'purple',srt = 45)
# abline(v=c(0), lwd=0.7, col="gray40", lty=2) # add central gridlines 
# abline(h=c(0), lwd=0.7, col="gray40",lty=2) # add central gridlines
# title("Changhsingian", line = 1)#, row=2
# Plot_ConvexHull(xcoord=PCA_data$x[changhsingian,][,1], ycoord= PCA_data$x[changhsingian,][,2], lcolor="#DB520F", bgcolor="#DB520F")


Plot_ConvexHull(xcoord=PCA_data$x[changhsingian,][,1], ycoord= PCA_data$x[changhsingian,][,2],
                lcolor="#FFA500", bgcolor="#FFA500","ALL",pch=0,D=TRUE)#FALSE

Plot_ConvexHull(xcoord=PCA_data$x[Induan,][,1], ycoord= PCA_data$x[Induan,][,2],
                lcolor="#0066FF", bgcolor="#0066FF","late Griesbachian-Dienerian",pch=2,D=TRUE)

# Plot_ConvexHull(xcoord=PCA_data$x[Lgriesbachian,][,1], ycoord= PCA_data$x[Lgriesbachian,][,2],
#                 lcolor="#EE82EE", bgcolor="#EE82EE","Lgriesbachian",pch=2,D=TRUE)
# Plot_ConvexHull(xcoord=PCA_data$x[Dienerian,][,1], ycoord= PCA_data$x[Dienerian,][,2],
#                 lcolor="#7B68EE", bgcolor="#7B68EE","Dienerian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[Smithian,][,1], ycoord= PCA_data$x[Smithian,][,2],
                lcolor="#F04028", bgcolor="#F04028","Smithian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[Spathian,][,1], ycoord= PCA_data$x[Spathian,][,2],
                lcolor="#F04028", bgcolor="#F04028","Smithian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[EAnisian,][,1], ycoord= PCA_data$x[EAnisian,][,2],
                lcolor="#0066FF", bgcolor="#0066FF","EAnisian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[LAnisian,][,1], ycoord= PCA_data$x[LAnisian,][,2], lcolor="#0066FF", bgcolor="#0066FF","LAnisian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[ELadinian,][,1], ycoord= PCA_data$x[ELadinian,][,2],
                lcolor="#9932CC", bgcolor="#9932CC","ELadinian",pch=2,D=TRUE)

Plot_ConvexHull(xcoord=PCA_data$x[LLadinian,][,1], ycoord= PCA_data$x[LLadinian,][,2],
                lcolor="#FF69B4", bgcolor="#FF69B4","LLadinian",pch=2,D=TRUE)


Time.groups <- list(changhsingian,LGriesbachian,Dienerian,Induan,Smithian,Spathian,EAnisian,LAnisian,ELadinian,LLadinian)
names(Time.groups) <- c("Changhsingian","Lgriesbachian","Dienerian","Induan","Smithian","Spathian","EAnisian","LAnisian","ELadinian","LLadinian")
#Time.groups
#Time.groups2=Time.groups[-1]##delete some columns
num_taxa_erp_interval = c(length(changhsingian),length(LGriesbachian),length(Dienerian),length(Induan),length(Smithian),length(Spathian),
                          length(EAnisian),length(LAnisian),length(ELadinian),length(LLadinian))
names(num_taxa_erp_interval) <- names(Time.groups) 
plot(num_taxa_erp_interval, type = "o", xlab="", ylab="num of lineages", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"),lwd=4,xaxt = "n")
axis(1, at = c(1:length(Time.groups)), labels = names(Time.groups))
text(num_taxa_erp_interval, labels = num_taxa_erp_interval, pos = 4, cex = 1)


################## different disparity metrics ##################
################### variances ###################
disparity.variances <- dispRity.per.group(PCA_data$x[,1:10], group = Time.groups, metric = c(sum, variances))
disparity.variances
summary(disparity.variances, digits = 5)
plot(disparity.variances, xlab="", ylab="Sum of variances D1-10", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))#, panel.first = grid(10,30),ylim=c(0,0.03)

hist(disparity.data$disparity$Changhsingian[2],col='skyblue',border = F)#breaks=28,
t_result <- t.test(disparity.data$disparity$Changhsingian[2], disparity.data$disparity$Induan[2], paired = TRUE)
t_result
################### distances ################### 这是对每个元素做的，因此需要做mean
disparity.centroids <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(mean, centroids))
disparity.centroids
summary(disparity.centroids, digits = 5)
plot(disparity.centroids , xlab="", ylab="mean of distances D1-2", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))

################### ranges ###################
disparity.ranges <- dispRity.per.group(PCA_data$x[,1:10], group = Time.groups, metric = c(sum, ranges))
disparity.ranges
summary(disparity.ranges, digits = 5)
plot(disparity.ranges , xlab="", ylab="Sum of ranges D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### diagonal ###################
disparity.diagonal <- dispRity.per.group(PCA_data$x[,1:10], group = Time.groups, metric = diagonal)
disparity.diagonal
summary(disparity.diagonal, digits = 5)
plot(disparity.diagonal , xlab="", ylab="Sum of diagonal D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### pairwise.dis ###################pairwise.dis
pairwise.disparity.dist <- dispRity.per.group(PCA_data$x[,1:10], group = Time.groups, metric = c(mean,pairwise.dist))
pairwise.disparity.dist
summary(pairwise.disparity.dist, digits = 5)
plot(pairwise.disparity.dist , xlab="", ylab="Sum of pairwise.dist D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### neighbours ###################
# neighbours <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(mean,neighbours))
# neighbours
# summary(neighbours, digits = 5)
# plot(neighbours , xlab="", ylab="Sum of neighbours D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### radius ###################
disparity.radius <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(sum,radius))
disparity.radius
summary(disparity.radius, digits = 5)
plot(disparity.radius , xlab="", ylab="Sum of radius D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### quantiles ###################
disparity.quantiles <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(sum,quantiles))
disparity.quantiles
summary(disparity.quantiles, digits = 5)
plot(disparity.quantiles , xlab="", ylab="Sum of quantiles D1-2", col=c( "#FF6666","#9966CC","#2bae85"))



################## order level validation##################Podocopida
# All orders: Ceratitida, Goniatitida, Phylloceratida, Prolecanitida
# Triassic: Ceratitida, Phylloceratida, Prolecanitida
pt_orders <- c('Ceratitida', 'Goniatitida', 'Phylloceratida', 'Prolecanitida')
triassic_orders <- c('Ceratitida', 'Phylloceratida', 'Prolecanitida')

Time.groups <- list(changhsingian,LGriesbachian,Dienerian,Induan,Smithian,Spathian,EAnisian,LAnisian,ELadinian,LLadinian)
names(Time.groups) <- c("Changhsingian","LGriesbachian","Dienerian","Induan","Smithian","Spathian","EAnisian","LAnisian","ELadinian","LLadinian")

Plot_ConvexHullsave<-function(xcoord, ycoord, lcolor, bgcolor,filename,pch,D,save = FALSE){
  par(new=D)
  plot(xcoord, ycoord,  xlab = "PC1", ylab = "PC2",pch = pch, cex=2,axes = !D, xlim = c(min(PCA_data$x[,1]),max(PCA_data$x[,1])), ylim = c(min(PCA_data$x[,2]),max(PCA_data$x[,2])))
  points(mean(xcoord),mean(ycoord), pch = pch+15, cex = 2, xlim = c(-0.3,0.5), ylim = c(min(PCA_data$x[,2]),max(PCA_data$x[,2])))
  cat(mean(xcoord),mean(ycoord))
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.4) , border = 0)
  title(filename)
}

filtered_rows_results <- list()

for (order_name in triassic_orders) {
  for (time_group in names(Time.groups)) {
    time_group <- as.character(time_group)
    
    filtered_rows <- which(Group$Order == order_name & Group[,time_group] == 1)
    
    #：Order_TimeGroup
    key_name <- paste(order_name, time_group, sep = "_")
    filtered_rows_results[[key_name]] <- filtered_rows
    
    if (length(filtered_rows) >2){
      Plot_ConvexHullsave(xcoord=PCA_data$x[changhsingian,][,1], ycoord= PCA_data$x[changhsingian,][,2],lcolor="#FFA500", bgcolor="#FFA500","",pch=0,D=FALSE,save=FALSE)#FALSE TRUE
      Plot_ConvexHullsave(xcoord=PCA_data$x[filtered_rows,][,1], ycoord= PCA_data$x[filtered_rows,][,2],lcolor="#0066FF", bgcolor="#0066FF",key_name,pch=2,D=TRUE,save=TRUE)
    }
    
    cat("\nOrder:", order_name, " Time Group:", time_group, 
        "\nFiltered Row Indices:\n", filtered_rows, "\n")
  }
}

dis.order.variances <- dispRity.per.group(PCA_data$x[,1:10], group = Filter(function(x) length(x) >= 3, filtered_rows_results), metric = c(sum, variances))
dis.order.variances
summary(dis.order.variances, digits = 5)
plot(dis.order.variances, xlab="", ylab="Sum of variances D1-10", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))#, panel.first = grid(10,30),ylim=c(0,0.03)

dis.order.ranges <- dispRity.per.group(PCA_data$x[,1:10], group = Filter(function(x) length(x) >= 3, filtered_rows_results), metric = c(sum, ranges))
dis.order.ranges
summary(dis.order.ranges, digits = 5)
plot(dis.order.ranges, xlab="", ylab="Sum of ranges D1-10", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))#, panel.first = grid(10,30),ylim=c(0,0.03)

lengths_info <- data.frame(Name = names(filtered_rows_results), 
                           Length = sapply(filtered_rows_results, length))
# save CSV
write.csv(lengths_info, "filtered_rows_lengths-ammonoids.csv", row.names = FALSE)



sum_of_variances=function(num,PCA_data,num_resampling){
  ran = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    ran[i] = sum(var(choice_pca[,1:10]))
    x_coord[i] = mean(xx)
    y_coord[i] = mean(yy)
  }
  return(ran)
}

sum_of_ranges=function(num,PCA_data,num_resampling){
  #a = n-m+1
  ran = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    ran[i] = sum(apply(choice_pca[,1:10], 2, function(x) max(x) - min(x)))
  }
  return(ran)
}

#######################################################
#################### PERMANOVA test####################
#######################################################

Interval_taxa<-function(bin){
  #taxa <- X.gpa$coords[,,bin]
  taxa <- PCA_data$x[bin,][,1:10]
  # COO2 <- aperm(taxa,c(2,1,3))
  # COO3 <- array(unlist(COO2),dim = c(length(bin),128))
  return (taxa)
}

changhsingian.taxa <- Interval_taxa(changhsingian)
PTTB.taxa <- Interval_taxa(PTTB)
Induan.taxa <- Interval_taxa(Induan)
Smithian.taxa <- Interval_taxa(Smithian)
Spathian.taxa <- Interval_taxa(Spathian)
EAnisian.taxa <- Interval_taxa(EAnisian)

total.taxa <- rbind(changhsingian.taxa, PTTB.taxa,Induan.taxa,Smithian.taxa,Spathian.taxa,EAnisian.taxa)
total.group<- c(rep("changhsingian",length(changhsingian)), rep("Induan",length(Induan)),rep("Smithian",length(Smithian)), rep("Spathian",length(Spathian)),rep("EAnisian",length(EAnisian)))

# total.taxa <- rbind(extinctions.taxa,newcomers.taxa)
# total.group<- c(rep("extinctions",length(extinctions)),rep("newcomers",length(newcomers)))

pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
  
{
  library(vegan)
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  for(elem in 1:ncol(co)){
    
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  
}

#pairwise.adonis(total.taxa, total.group, sim.method="euclidean", p.adjust.m= "bonferroni")

total.group2 <- c(rep("changhsingian",length(changhsingian)),rep("EAnisian",length(EAnisian)))
xx = rbind(changhsingian.taxa[,1:10],EAnisian.taxa[,1:10])

dis <- vegdist(xx,method = 'euclidean')# gower euclidean canberra
mod = betadisper(dis,total.group2)#, type = c("median","centroid")
permutest(mod)

#######################################################
################## random simulations##################
#######################################################
range=function(num,PCA_data,num_resampling){
  #a = n-m+1
  ran = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    ran[i] = sum(apply(choice_pca[,1:10], 2, function(x) max(x) - min(x)))
    # xx = choice_pca[,1]
    # yy = choice_pca[,2]
    # #yy = sample(PTA$Comp2,size=num,replace=TRUE)
    # ran[i] = (max(xx)-min(xx))+(max(yy)-min(yy))
    #ran[i] = var(xx)+var(yy)
  }
  return(ran)
}

variance=function(num,PCA_data,num_resampling){
  #a = n-m+1
  ran = rep(NA,num_resampling)
  x_coord = rep(NA,num_resampling)
  y_coord = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    ran[i] = sum(apply(choice_pca[,1:10], 2, var))
  return(ran)
}

SQS_div <- c(21,23,90,69,49)
times <- list(changhsingian, Induan,Olenekian,	Anisian,	Ladinian)
xaxes  <- c('changhsingian', 'Induan','Olenekian',	'Anisian',	'Ladinian')
num_resampling <- 1000
total_variance <- matrix(NA, nrow = 5, ncol = 1000)
total_range <- matrix(NA, nrow = 5, ncol = 1000)
for (i in 1:5){
  num <- SQS_div[i]
  time <- times[i]
  total_variance[i,] <- variance(num,PCA_data$x[times[[i]],],num_resampling)
  total_range[i,] <- range(num,PCA_data$x[times[[i]],],num_resampling)
}

output_variances = c(apply(total_variance, 1, function(x) quantile(x, 0.025)),apply(total_variance, 1, mean),apply(total_variance, 1, function(x) quantile(x, 0.975)))
output_variances

output_ranges = c(apply(total_range, 1, function(x) quantile(x, 0.025)),apply(total_range, 1, mean),apply(total_range, 1, function(x) quantile(x, 0.975)))
output_ranges

