library(geomorph)
library(dispRity)
library(factoextra)
library(ggpubr)

df = read.csv("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/ammonoid/ammonoid_a_landmarks.csv",header=F,row.names = 1)
ammonoid.links <- as.matrix(read.table("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/ammonoid/ammonoid_links.txt"))
curves<-as.matrix(read.csv("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/ammonoid/ammonoid_curveslide.csv",header = T))
num_sp = length(row.names(df))
AA <- array(data = unlist(t(df[,1:128])),dim = c(2,64,num_sp),dimnames=NULL)
X.gpa <- gpagen(BB, curves=curves,PrinAxes = FALSE,Proj=T)
# summary(X.gpa)
# plot(X.gpa)

average.shape <- mshape(X.gpa$coords)#average shape points
PCA_data<-gm.prcomp(X.gpa$coords)
plot(PCA_data)
summary(PCA_data)
plotRefToTarget(PCA_data$shapes$shapes.comp1$max,average.shape)#, links = ostracod.links,label=FALSE
plotAllSpecimens(X.gpa$coords,links = ammonoid.links)#

pvar <- (PCA_data$sdev^2)/(sum(PCA_data$sdev^2))*100
names(pvar) <- seq(1:length(pvar))
par(mgp=c(2,1,0))
barplot(pvar[1:10], col=rgb(0.2,0.4,0.6,0.6),ylim=c(0,80), xlab= "Principal Components", ylab = "% Variance") # save this plot
title("Proportion of variance", line = 1)#, row=2
pvar[1:3]

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.27) , border = 0)
}

Group <- read.csv("C:/Users/xiaokang/Jupyter_code/feature_extraction/PTimages/landmarks/ammonoid/Ammonoidea_interval.csv",row.names = 1)
families <- Group[,'Order']

changhsingian <- row.names(Group)[Group[,'Changhsingian'] == "1"]
PTTB <- row.names(Group)[Group[,'PTTB'] == "1"]
Lgriesbachian <- row.names(Group)[Group[,'LGriesbachian'] == "1"]
#Griesbachian <- row.names(Group)[Group[,'Griesbachian'] == "1"]
Dienerian <- row.names(Group)[Group[,'Dienerian'] == "1"]
lGriDie <- row.names(Group)[Group[,'lGriDie'] == "1"]
#Induan <- row.names(Group)[Group[,'Induan'] == "1"]
survivors <- row.names(Group)[Group[,'survivors'] == "1"]
extinctions <- row.names(Group)[Group[,'extinctions'] == "1"]
newcomers <- row.names(Group)[Group[,'newcomers'] == "1"]

plot(PCA_data$x[,1], PCA_data$x[,2],  xlab = "PC1", ylab = "PC2", pch = 20, cex=1)
text(x = PCA_data$x[,1], y = PCA_data$x[,2],labels = seq(1:num_sp),cex = 1,col = 'purple',srt = 45)
abline(v=c(0), lwd=0.7, col="gray40", lty=2) # add central gridlines 
abline(h=c(0), lwd=0.7, col="gray40",lty=2) # add central gridlines
#title("all species")

Plot_ConvexHull<-function(xcoord, ycoord, lcolor, bgcolor,interval, pch,D){
  par(new=D)
  #plot(xcoord, ycoord,  xlab = "PC1", ylab = "PC2",pch = pch, cex=1,axes = !D, xlim = c(-0.4,0.25), ylim = c(-0.15,0.1))
  points(mean(xcoord),mean(ycoord), pch = pch+15, cex = 2)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)
  cat(mean(xcoord),mean(ycoord))
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(x = xcoord[hpts], y =  ycoord[hpts], col = adjustcolor(bgcolor, alpha.f = 0.4) , border = 0)
  #title(interval)
  
}   

#changhsingian
Plot_ConvexHull(xcoord=PCA_data$x[changhsingian,][,1], ycoord= PCA_data$x[changhsingian,][,2],
                lcolor="#F04028", bgcolor="#F04028","Changhsingian",pch=0,D=FALSE)
#points(mean(PCA_data$x[changhsingian,][,1]),mean(PCA_data$x[changhsingian,][,2]), xlim = c(-0.42,0.3), ylim = c(-0.35,0.2), pch = 15, cex = 2)
Plot_ConvexHull(xcoord=PCA_data$x[PTTB,][,1], ycoord= PCA_data$x[PTTB,][,2],
                lcolor="#696969", bgcolor="#696969","early Griesbachian",pch=1,D=TRUE)
#points(mean(PCA_data$x[PTTB,][,1]),mean(PCA_data$x[PTTB,][,2]), xlim = c(-0.42,0.3), ylim = c(-0.35,0.2), pch = 20, cex = 2)
# Plot_ConvexHull(xcoord=PCA_data$x[Lgriesbachian,][,1], ycoord= PCA_data$x[Lgriesbachian,][,2],
#                 lcolor="#0066FF", bgcolor="#0066FF","late Griesbachian",pch=2,D=TRUE)
# #Plot_ConvexHull(xcoord=PCA_data$x[Griesbachian,][,1], ycoord= PCA_data$x[Griesbachian,][,2], lcolor="#0066FF", bgcolor="#0066FF","Griesbachian",pch=2)
# Plot_ConvexHull(xcoord=PCA_data$x[Dienerian,][,1], ycoord= PCA_data$x[Dienerian,][,2], 
#                 lcolor="#000000", bgcolor="#FFFF00","Dienerian",pch=5,D=TRUE)
Plot_ConvexHull(xcoord=PCA_data$x[lGriDie,][,1], ycoord= PCA_data$x[lGriDie,][,2],
                lcolor="#0066FF", bgcolor="#0066FF","late Griesbachian-Dienerian",pch=2,D=TRUE)

points(mean(PCA_data$x[survivors,][,1]),mean(PCA_data$x[survivors,][,2]), pch = 3, cex = 2)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.2)
points(mean(PCA_data$x[newcomers,][,1]),mean(PCA_data$x[newcomers,][,2]), pch = 4, cex = 2)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.2)
cat(mean(PCA_data$x[survivors,][,1]),mean(PCA_data$x[survivors,][,2]))
cat(mean(PCA_data$x[newcomers,][,1]),mean(PCA_data$x[newcomers,][,2]))
# Plot_ConvexHull(xcoord=PCA_data$x[survivors,][,1], ycoord= PCA_data$x[survivors,][,2], lcolor="#AC5500BB", bgcolor="#AC5500BB","Survivors",pch=5)
# Plot_ConvexHull(xcoord=PCA_data$x[newcomers,][,1], ycoord= PCA_data$x[newcomers,][,2], lcolor="#F3D32C", bgcolor="#F3D32C","newcomers",pch=1)
# Plot_ConvexHull(xcoord=PCA_data$x[extinctions,][,1], ycoord= PCA_data$x[extinctions,][,2], lcolor="#444444", bgcolor="#444444","extinctions",pch=0)

var_cha = var(PCA_data$x[changhsingian,][,1])+var(PCA_data$x[changhsingian,][,2])#+var(PCA_data$x[changhsingian,][,3])
#var_gri = var(PCA_data$x[Griesbachian,][,1])+var(PCA_data$x[Griesbachian,][,2])#+var(PCA_data$x[Griesbachian,][,3])
var_egri = var(PCA_data$x[PTTB,][,1])+var(PCA_data$x[PTTB,][,2])#+var(PCA_data$x[Egriesbachian,][,3])
var_lgri = var(PCA_data$x[Lgriesbachian,][,1])+var(PCA_data$x[Lgriesbachian,][,2])#+var(PCA_data$x[Lgriesbachian,][,3])
var_die = var(PCA_data$x[Dienerian,][,1])+var(PCA_data$x[Dienerian,][,2])#+var(PCA_data$x[Dienerian,][,3])
var_lGriDie = var(PCA_data$x[lGriDie,][,1])+var(PCA_data$x[lGriDie,][,2])#+var(PCA_data$x[lGriDie,][,3])
var_cha 
#var_gri
var_egri
var_lgri
var_die
var_lGriDie


################## disparity metrics ##################
# disparity calculation
####Disparity
#Time.groups <- list(changhsingian,Egriesbachian,Lgriesbachian,Dienerian,newcomers)# volume is inapplicable for survivors 
#names(Time.groups) <- c("Changhsingian","Egriesbachian","Lgriesbachian","Dienerian","newcomers")
#Time.groups
Time.groups <- list(changhsingian,PTTB,Lgriesbachian,Dienerian,lGriDie,survivors,extinctions,newcomers)
names(Time.groups) <- c("Changhsingian","PTTB","Lgriesbachian","Dienerian","lGriDie","survivors","extinctions","newcomers")
#Time.groups
#Time.groups2=Time.groups[-1]##delete some columns
num_taxa_erp_interval = c(length(changhsingian),length(PTTB),length(Lgriesbachian),length(Dienerian),length(lGriDie),length(survivors),length(extinctions),length(newcomers))
names(num_taxa_erp_interval) <- names(Time.groups) 
plot(num_taxa_erp_interval, type = "o", xlab="", ylab="num of lineages", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"),lwd=4)

################### variances ###################
disparity.data <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(sum, variances))
disparity.data
summary(disparity.data, digits = 5)
plot(disparity.data, xlab="", ylab="Sum of variances D1-2", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))#, panel.first = grid(10,30),ylim=c(0,0.03)

hist(disparity.data$disparity$Changhsingian[2],col='skyblue',border = F)#breaks=28,
t_result <- t.test(disparity.data$disparity$Changhsingian[2], disparity.data$disparity$lGriDie[2], paired = TRUE)
t_result
################### distances ################### 这是对每个元素做的，因此需要做mean
disparity.centroids <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(mean, centroids))
disparity.centroids
summary(disparity.centroids, digits = 5)
plot(disparity.centroids , xlab="", ylab="mean of distances D1-2", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))

################### ranges ###################
ranges.disparity.data <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(sum, ranges))
ranges.disparity.data
summary(ranges.disparity.data, digits = 5)
plot(ranges.disparity.data , xlab="", ylab="Sum of ranges D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### diagonal ###################
diagonal.disparity.data <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = diagonal)
diagonal.disparity.data
summary(diagonal.disparity.data, digits = 5)
plot(diagonal.disparity.data , xlab="", ylab="Sum of diagonal D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

################### pairwise.dis ###################pairwise.dis
pairwise.disparity.dist <- dispRity.per.group(PCA_data$x[,1:2], group = Time.groups, metric = c(mean,pairwise.dist))
pairwise.disparity.dist
summary(pairwise.disparity.dist, digits = 5)
plot(pairwise.disparity.dist , xlab="", ylab="Sum of pairwise.dist D1-2", col=c( "#FF6666","#9966CC","#2bae85"))

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


################## density plots of survivors, extinctions, and newcomers##################
classes = rep(NA,num_sp)
num <- 1:num_sp
for(i in num){
  if (Group[,'survivors'][i]==1){
    classes[i] = "survivor"
  }
  else if(Group[,'extinctions'][i]==1){
    classes[i] = "extinction"
  }
  else if(Group[,'newcomers'][i]==1){
    classes[i] = "newcomer"
  }
}
#classes
data2 <- data.frame(x=PCA_data$x[,1],y=PCA_data$x[,2],class=classes)

ggscatterhist(
  data2,  x ='x', y = 'y',  xlab = ' PC 1',ylab = 'PC 2',
  shape=21,color ="black",fill= "class", size =3, alpha = 0.8,
  palette = c("#FC4E07", "#E7B800", "#00AFBB"),
  margin.plot =  "density",#histogram  density  boxplot
  margin.params = list(fill = "class", color = "black", size = 0.2),
  legend = c(0.9,0.15),
  ggtheme = theme_minimal())


################## random extinction simulations ##################
variance=function(num,PCA_data,num_resampling){
  ran = rep(NA,num_resampling)
  x_coord = rep(NA,num_resampling)
  y_coord = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    xx = choice_pca[,1]
    yy = choice_pca[,2]
    #yy = sample(PTA$Comp2,size=num,replace=TRUE)
    #ran[i] = (max(xx)-min(xx))+(max(yy)-min(yy))
    ran[i] = var(xx)+var(yy)
    x_coord[i] = mean(xx)
    y_coord[i] = mean(yy)
  }
  # cat(mean(x_coord),quantile(x_coord,0.025),quantile(x_coord,0.975))
  # cat(mean(y_coord),quantile(y_coord,0.025),quantile(y_coord,0.975))
  return(ran)
}

range=function(num,PCA_data,num_resampling){
  ran = rep(NA,num_resampling)
  for (i in 1:num_resampling){
    
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    xx = choice_pca[,1]
    yy = choice_pca[,2]
    ran[i] = (max(xx)-min(xx))+(max(yy)-min(yy))
    #ran[i] = var(xx)+var(yy)
  }
  return(ran)
}

###################### resampling num_survivors ######################
survivors_ranges <- (max(PCA_data$x[survivors,][,1])-min(PCA_data$x[survivors,][,1]))+(max(PCA_data$x[survivors,][,2])-min(PCA_data$x[survivors,][,2]))
survivors_ranges
survivors_variances <- var(PCA_data$x[survivors,][,1])+var(PCA_data$x[survivors,][,2])
survivors_variances

num_survivors = length(survivors)
num_changhsingian = length(changhsingian)
num_resampling = 10000
sub_changhsingian = variance(num_survivors,PCA_data$x[changhsingian,],num_resampling)
output_cha = c(mean(sub_changhsingian),sd(sub_changhsingian),quantile(sub_changhsingian,0.025),quantile(sub_changhsingian,0.975))
output_cha

sub_survivors = variance(num_survivors,PCA_data$x[survivors,],num_resampling)
output_sur = c(mean(sub_survivors),sd(sub_survivors),quantile(sub_survivors,0.025),quantile(sub_survivors,0.975))
hist(sub_changhsingian,col='skyblue',border = F,breaks=56,xlim = c(0,0.035),ylim = c(0,800))#breaks=28,
#hist(sub_survivors, col='#00AFBB',add=T,breaks=28,border = F)#, col=rgb(1,0,0,1/4), xlim=c(0,0.08),breaks=28
#abline(v=c(output_cha[1]), lwd=1, col="gray10",lty=1) # add central gridlines
abline(v=c(survivors_variances), lwd=3, col="red",lty=1) # add central gridlines
abline(v=c(quantile(sub_changhsingian,0.05)), lwd=3, col="gray40",lty=1)

###################### resampling num_lGriDie ######################
lGriDie_variances <- var(PCA_data$x[lGriDie,][,1])+var(PCA_data$x[lGriDie,][,2])
lGriDie_variances

num_lGriDie=length(lGriDie)
sub_changhsingian2 = variance(num_lGriDie,PCA_data$x[changhsingian,],num_resampling)
output_cha2 = c(mean(sub_changhsingian2),sd(sub_changhsingian2),quantile(sub_changhsingian2,0.025),quantile(sub_changhsingian2,0.25),quantile(sub_changhsingian2,0.75),quantile(sub_changhsingian2,0.975))
output_cha2
sub_griesbachian = variance(num_lGriDie,PCA_data$x[lGriDie,],num_resampling)
output_gri = c(mean(sub_griesbachian),sd(sub_griesbachian),quantile(sub_griesbachian,0.025),quantile(sub_changhsingian2,0.25),quantile(sub_changhsingian2,0.75),quantile(sub_griesbachian,0.975))
output_gri
hist(sub_changhsingian2, col="skyblue",breaks=56,border = F,xlim = c(0,0.04),ylim = c(0,800))#col=rgb(0,0,1,1/4)
abline(v=c(lGriDie_variances), lwd=3, col="red",lty=1) # add central gridlines
abline(v=c(quantile(sub_changhsingian2,0.05)), lwd=3, col="gray40",lty=1)

sub_changhsingian2 = range(num_lGriDie,PCA_data$x[changhsingian,],num_resampling)
output_cha2 = c(mean(sub_changhsingian2),sd(sub_changhsingian2),quantile(sub_changhsingian2,0.025),quantile(sub_changhsingian2,0.25),quantile(sub_changhsingian2,0.75),quantile(sub_changhsingian2,0.975))
output_cha2

center=function(num,PCA_data){
  #print(PCA_data[,1])
  center_x = mean(PCA_data[,1])
  center_y = mean(PCA_data[,2])
  distance = rep(NA,10000)
  for (i in 1:10000){
    
    choice = sample(row.names(PCA_data),size=num,replace=TRUE)
    choice_pca = PCA_data[choice,]
    xx = mean(choice_pca[,1])
    yy = mean(choice_pca[,2])
    distance[i] = abs(center_x-xx)+abs(center_y-yy)#((center_x-xx)^2+(center_y-yy)^2)^0.5 #
  }
  return(distance)
}

num_lGriDie = length(lGriDie)
random_extinction = center(num_lGriDie,PCA_data$x[changhsingian,])
# observed_dis = ((mean(PCA_data$x[changhsingian,][,1])-mean(PCA_data$x[lGriDie,][,1]))^2+(mean(PCA_data$x[changhsingian,][,2])-mean(PCA_data$x[lGriDie,][,2]))^2)^0.5
# observed_dis
observed_cent_rang = abs(mean(PCA_data$x[changhsingian,][,1])-mean(PCA_data$x[lGriDie,][,1]))+abs(mean(PCA_data$x[changhsingian,][,2])-mean(PCA_data$x[lGriDie,][,2]))
observed_cent_rang
output_random_extinction = c(mean(random_extinction),sd(random_extinction),quantile(random_extinction,0.025),quantile(random_extinction,0.975))
output_random_extinction
hist(random_extinction, col=rgb(0,0,1,1/4),breaks=28,border = F,xlim = c(0,0.17))
abline(v=c(quantile(random_extinction,0.95)), lwd=3, col="gray40",lty=2)
abline(v=observed_dis, lwd=3, col="red",lty=2)


################################ PERMANOVA test################################
Interval_taxa<-function(bin){
  taxa <- PCA_data$x[bin,][,1:3]
  return (taxa)
}

changhsingian.taxa <- Interval_taxa(changhsingian)
PTTB.taxa <- Interval_taxa(PTTB)
lGriDie.taxa <- Interval_taxa(lGriDie)
extinctions.taxa <- Interval_taxa(extinctions)
newcomers.taxa <- Interval_taxa(newcomers)
survivors.taxa <- Interval_taxa(survivors)

total.taxa <- rbind(changhsingian.taxa, PTTB.taxa,lGriDie.taxa,extinctions.taxa,newcomers.taxa)
total.group<- c(rep("changhsingian",length(changhsingian)), rep("PTTB",length(PTTB)),rep("lGriDie",length(lGriDie)), rep("extinctions",length(extinctions)),rep("newcomers",length(newcomers)))

total.group2 <- c(rep("changhsingian",lengthchanghsingian),rep("lGriDie",length(lGriDie)))
total_rbind = rbind(changhsingian.taxa[,1:3],lGriDie.taxa[,1:3])
dis <- vegdist(total_rbind,method = 'euclidean')# gower euclidean canberra
mod = betadisper(dis,total.group2)#, type = c("median","centroid")
permutest(mod)


# pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
# {
#   library(vegan)
#   co = as.matrix(combn(unique(factors),2))
#   pairs = c()
#   F.Model =c()
#   R2 = c()
#   p.value = c()
#   for(elem in 1:ncol(co)){
#     ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
#                   factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
#     pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
#     F.Model =c(F.Model,ad$aov.tab[1,4]);
#     R2 = c(R2,ad$aov.tab[1,5]);
#     p.value = c(p.value,ad$aov.tab[1,6])
#   }
#   p.adjusted =p.adjust(p.value,method=p.adjust.m)
#   pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
#   return(pairw.res)
# }
# pairwise.adonis(total.taxa, total.group, sim.method="euclidean", p.adjust.m= "bonferroni")


################## family level validation##################
# plot(PCA_data$x[,1], PCA_data$x[,2],  xlab = "PC1", ylab = "PC2", pch = 20, cex=1)
# text(x = PCA_data$x[,1], y = PCA_data$x[,2],labels = seq(1:num_sp),cex = 1,col = 'purple',srt = 45)
# abline(v=c(0), lwd=0.7, col="gray40", lty=2) # add central gridlines 
# abline(h=c(0), lwd=0.7, col="gray40",lty=2) # add central gridlines
# #title("all species")
# 
# Xenodiscidae <- intersect(changhsingian,row.names(Group)[Group[,'Family'] == "Xenodiscidae"])
# Dzhulfitidae <- intersect(changhsingian,row.names(Group)[Group[,'Family'] == "Dzhulfitidae"])
# Araxoceratidae <- intersect(changhsingian,row.names(Group)[Group[,'Family'] == "Araxoceratidae"])
# Pleuronodoceratidae <- intersect(changhsingian,row.names(Group)[Group[,'Family'] == "Pleuronodoceratidae"])
# Medlicottiidae <- row.names(Group)[Group[,'Family'] == "Medlicottiidae"]
# 
# Liuchengoceratidae <- row.names(Group)[Group[,'Family'] == "Liuchengoceratidae"]
# Paragastrioceratidae <- row.names(Group)[Group[,'Family'] == "Paragastrioceratidae"]
# Pleuronodoceratidae <- row.names(Group)[Group[,'Family'] == "Pleuronodoceratidae"]
# Tapashanitidae <- row.names(Group)[Group[,'Family'] == "Tapashanitidae"]
# 
# #Families
# Plot_ConvexHull(xcoord=PCA_data$x[Xenodiscidae,][,1], ycoord= PCA_data$x[Xenodiscidae,][,2],
#                 lcolor="#DB520F", bgcolor="#DB520F","Changhsingian",pch=0)##pink
# Plot_ConvexHull(xcoord=PCA_data$x[Dzhulfitidae,][,1], ycoord= PCA_data$x[Dzhulfitidae,][,2],
#                 lcolor="#444444", bgcolor="#444444","early Griesbachian",pch=2)##grey
# Plot_ConvexHull(xcoord=PCA_data$x[Araxoceratidae,][,1], ycoord= PCA_data$x[Araxoceratidae,][,2],
#                 lcolor="#0066FF", bgcolor="#0066FF","late Griesbachian",pch=2)#blue
# Plot_ConvexHull(xcoord=PCA_data$x[Pleuronodoceratidae,][,1], ycoord= PCA_data$x[Pleuronodoceratidae,][,2], 
#                 lcolor="#F3D32C", bgcolor="#F3D32C","Dienerian",pch=1)##yellow

# Plot_ConvexHull(xcoord=PCA_data$x[Liuchengoceratidae,][,1], ycoord= PCA_data$x[Liuchengoceratidae,][,2],
#                 lcolor="#DB520F", bgcolor="#DB520F","Changhsingian",pch=0)##pink
# Plot_ConvexHull(xcoord=PCA_data$x[Paragastrioceratidae,][,1], ycoord= PCA_data$x[Paragastrioceratidae,][,2],
#                 lcolor="#444444", bgcolor="#444444","early Griesbachian",pch=2)##grey
# Plot_ConvexHull(xcoord=PCA_data$x[Pleuronodoceratidae,][,1], ycoord= PCA_data$x[Pleuronodoceratidae,][,2],
#                 lcolor="#0066FF", bgcolor="#0066FF","late Griesbachian",pch=2)#blue
# Plot_ConvexHull(xcoord=PCA_data$x[Tapashanitidae,][,1], ycoord= PCA_data$x[Tapashanitidae,][,2], 
#                 lcolor="#F3D32C", bgcolor="#F3D32C","Dienerian",pch=1)##yellow

# points(mean(PCA_data$x[Liuchengoceratidae,][,1]),mean(PCA_data$x[Liuchengoceratidae,][,2]), pch = 0+15, cex = 2,col='grey')
# points(mean(PCA_data$x[Paragastrioceratidae,][,1]),mean(PCA_data$x[Paragastrioceratidae,][,2]), pch = 0+15, cex = 2,col='grey')
# points(mean(PCA_data$x[Pleuronodoceratidae,][,1]),mean(PCA_data$x[Pleuronodoceratidae,][,2]), pch = 0+15, cex = 2,col='grey')
# points(mean(PCA_data$x[Tapashanitidae,][,1]),mean(PCA_data$x[Tapashanitidae,][,2]), pch = 0+15, cex = 2,col='grey')
# 
# points(mean(PCA_data$x[Medlicottiidae,][,1]),mean(PCA_data$x[Medlicottiidae,][,2]), pch = 0+15, cex = 2,col='red')#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)
# points(mean(PCA_data$x[Otoceratidae,][,1]),mean(PCA_data$x[Otoceratidae,][,2]), pch = 2+15, cex = 2,col='red')#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)
# points(mean(PCA_data$x[changhsingian,][,1]),mean(PCA_data$x[changhsingian,][,2]), pch = 2+15, cex = 4)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)
# points(mean(PCA_data$x[survivors,][,1]),mean(PCA_data$x[survivors,][,2]), pch = 1+15, cex = 3)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)
# points(mean(PCA_data$x[Griesbachian,][,1]),mean(PCA_data$x[Griesbachian,][,2]), pch = 2+15, cex = 4)#, xlim = c(-0.42,0.3), ylim = c(-0.35,0.3)

# Family.groups <- list(Xenodiscidae,Dzhulfitidae,Araxoceratidae,Pleuronodoceratidae)
# names(Family.groups) <- c("Xenodiscidae","Dzhulfitidae","Araxoceratidae","Pleuronodoceratidae")
# 
# disparity.family <- dispRity.per.group(PCA_data$x[,1:2], group = Family.groups, metric = c(sum, variances))
# disparity.family
# summary(disparity.family, digits = 5)
# plot(disparity.family , xlab="", ylab="Sum of variances D1-2", col=c( "#FF6666","#9966CC","#2bae85","#FF8C00"))
