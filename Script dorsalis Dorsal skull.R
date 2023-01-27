########################################################################
########################################################################
#R script of "Geographic variation in select species of the bat genus Platyrrhinus"
#by Velazco, Ly, McAllister and Esquivel
# THERYA, 2023, Vol 14(1):121-130. DOI:10.12933/therya-23-2208

# This Script contains the following functions:
## Load raw data coordinates and classifiers
## GPA
## PCA and LDA/CVA
## Allometry 
## Size and shape analysis
## Graphics

#########################################################################
########################################################################
# Load packages
# install.packages("")
dir()
require(geomorph)
require(MASS)
require(ggplot2)
require(car)
require(vegan)
require(Morpho)

########################################################################
# Upload data (.tps)  
load("dorsalis_Dorsal.RData")

dim(tps)
dim(plan)
groups<-as.factor(plan[,2])
sex<-as.factor(plan[,3])

########################################################################
# GPA (Generalized Procrustes Analysis)
gpa.object<-gpagen(tps)
shape<-gpa.object$coords 
size<-gpa.object$Csize

plotAllSpecimens(shape)

########################################################################
# Find outliers. This step can be skipped
plotOutliers(shape) 
 plotOutliers(shape,groups = groups,inspect.outliers=TRUE)

# Remove specimens/outliers.If you want to remove one or more individuals from your dataset: here we want to remove individuals 127 and 131. 
tps<-tps[,,-c(127,131)] 
plan<-plan[-c(127,131),]

########################################################################

# Shape visualization
# Mean shape
ref<-mshape(shape)

# Define links between landmarks 
# links<-define.links(ref,ptsize=1)# Manually
links<-
  matrix(c(1,3,3,4,4,5,5,6,6,7,1,2),nrow=6,ncol=2,byrow=T)

# Plot specimens against mean shape
GP1<-gridPar(pt.bg="gray",link.col="gray",link.lty=1) 
plotRefToTarget(ref,shape[,,11],links=links,method="TPS",mag = 1) 
plotRefToTarget(ref,shape[,,11],links=links,method="vector", mag = 1) 
plotRefToTarget(ref,shape[,,11],links=links,method="points",gridPars=GP1) # target = black, reference = gray

plotAllSpecimens(shape,mean=TRUE,links=links)

############################ Size Analyses ##############################

plan$CS<-size # size per each group
chocoensis<-plan[plan$Species=="chocoensis",]
dorsalis<-plan[plan$Species=="dorsalis",]

##### 1.1 Does the size differ between Females and Males? (considering all as one group) ####
boxplot(size~sex,ylab="CENTROID SIZE")
#Vizualization of size differences between sexes
boxplot(log(size)~sex,ylab="LOG CENTROID SIZE")#quienes son mas grandes? macho o hembras? esto solo es para visualizar, aun no se ha hecho ningun test
#tiff("1. Boxplot log (size) vs sex.tiff", units="in", width=5, height=5, res=300)
#dev.off()

# Test
shapiro.test(log(size))
leveneTest(log(size),group = sex,center = "median")
t.test(log(size)~sex,var.equal = TRUE)

##### 1.2 Does the size differ between Females and Males inside each subpopulation or groups/species? (2 groups) ####
par(mfrow=c(1,2))
boxplot(log(chocoensis$CS)~chocoensis$Sex,ylab="LOG CENTROID SIZE")
boxplot(log(dorsalis$CS)~dorsalis$Sex,ylab="LOG CENTROID SIZE")

# Test: Variation in size inside species (considering two separate group)
shapiro.test(log(chocoensis$CS))
leveneTest(log(chocoensis$CS),group = chocoensis$Sex,center = "median")
t.test(log(chocoensis$CS)~chocoensis$Sex,var.equal = TRUE)

shapiro.test(log(dorsalis$CS))
leveneTest(log(dorsalis$CS),group = dorsalis$Sex,center = "median")
t.test(log(dorsalis$CS)~dorsalis$Sex,var.equal = TRUE)

##### 1.3 Sexual dimorphism in size? ####
boxplot(size~groups*sex,ylab="CENTROID SIZE")
boxplot(log(size)~groups*sex,ylab="LOG CENTROID SIZE")

# Test: ANOVA with interaction
# Unbalanced designs?
anova<-aov(log(size)~sex*groups)
Anova(anova, type = "II")
TukeyHSD(anova)

gdf <- geomorph.data.frame(shape= shape,size=size,
                           sex=plan$Sex,species=groups)

fit.size.sp.int <- lm.rrpp(size ~ sex * species,iter = 10000,RRPP = TRUE,
                           SS.type = "II", data = gdf, 
                           print.progress = FALSE) # Does size di???er between species, while accounting for size covarying with sexes?

anova(fit.size.sp.int)

boxplot(log(size)~groups,ylab="LOG CENTROID SIZE")
#tiff("Boxplot log (size) vs species.tiff", units="in", width=5, height=5, res=300)
#dev.off()

######################### Shape Analyses ################################

#### 2.1 Shape variation between species and sexes: sexual dimorphism? ####
new_shape<-two.d.array(shape)
gdf1 <- geomorph.data.frame(shape= new_shape,size=size,
                            sex=plan$Sex,species=groups)
fit.full.model <-lm.rrpp(shape ~ size * sex * species,iter = 10000,RRPP = TRUE,
                         SS.type = "II", data = gdf1, 
                         print.progress = FALSE)

anova(fit.full.model,effect.type = "F")

#### PCA-Normal-WITHOUT SEXUAL DIMORPHISM #### 
### Customize PCA 
col.group<-c("#010305","#4271AE") 
names(col.group)<-levels(groups)
col.group<-col.group[match(groups,names(col.group))]

PCA<-gm.prcomp(shape) 
xlab<-"Principal Component 1 (48.07%)"
ylab<-"Principal Component 2 (16.38%)"
mat<-matrix(c(4,5,0,1,1,2,1,1,3),3) # Split the plot window
layout(mat, widths=c(3,2,2), heights=c(1,1,1))
par(mar=c(4, 4, 1, 1))

plot(PCA$x[,1],PCA$x[,2],pch=21,cex=3.5,cex.lab=1.8,bg=col.group,xlab=xlab,ylab=ylab,font.lab = 2,font.axis = 2,asp=T)
plotRefToTarget(ref,PCA$shapes$shapes.comp1$min,links=links,method="points",gridPars=GP1,mag = 2)
plotRefToTarget(ref,PCA$shapes$shapes.comp1$max,links=links,method="points",gridPars=GP1,mag = 2)
plotRefToTarget(ref,PCA$shapes$shapes.comp2$max,links=links,method="points",gridPars=GP1,mag = 2)
plotRefToTarget(ref,PCA$shapes$shapes.comp2$min,links=links,method="points",gridPars=GP1,mag = 2)

par(mfrow=c(1,1))   

# 95% confidence ellipses
plot(PCA$x[,1],PCA$x[,2],pch=21,cex=2,bg=col.group,xlab=xlab,ylab=ylab,font.lab = 2,asp=T)
ordihull(PCA$x,group=groups,lwd = 1.5)

##### LDA / CVA ####

# LDA (Linear Discriminant Analysis) 
cva<-lda(PCA$x[,1:6],groups)
plot(cva)
cva<-lda(PCA$x[,1:6],groups,CV=T) #LDA com Jackknife cross validation
tab<-table(groups,cva$class) 
lda.p<-diag(tab)/summary(groups)*100
lda.p # providing correct classification for each group

##### Customize Figures #########
require(ggplot2)
fill <- c("#56B4E9","#4271AE")
Box_1<-ggplot(plan, aes(x=plan$Sex, y=log(plan$CS))) + 
  geom_boxplot(fill = fill,size = 1 )+
  labs(title="Dorsal view: plot of Centroid Size by sex",x="Sex", y = "Log (Centroid Size)")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold", hjust = 0.5),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9))

ggsave("1. Boxplot log (size) vs sex.tiff",units="in",width=5, height=5, dpi=300)
#
Box_2<-ggplot(plan, aes(x=groups, y=log(size))) + 
  geom_boxplot(fill = fill,size = 1 )+
  labs(title=NULL,x="Species", y = "Log (Centroid Size)")+
  theme_classic()+
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold", hjust = 0.5),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(colour="black", size = 11),
        axis.text.y = element_text(colour="black", size = 9))
Box_2 + scale_x_discrete(labels=c("Platyrrhinus chocoensis", "Platyrrhinus dorsalis"))

ggsave("2. Boxplot log (size) vs species.tiff",units="in",width=5, height=5, dpi=300)