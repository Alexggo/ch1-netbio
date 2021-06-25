
#### With networkdist/k.edge.connectivity
connect <- read.csv(file.path("data/5.Targets_NetworkDistance","EcoliNet.v1.benchmark.txtnetvalues_target.csv"))
connect <- connect %>% mutate(drugdrug=paste0(drug1,"-",drug2))



mat1 <- full_join(mat,connect,by="drugdrug")
mat2 <- mat1 %>% 
  filter(!is.na(ebw)&!is.na(ecr)&!is.na(seo)&!is.na(stm)&!is.na(pae)&!is.na(pau)) %>% 
  filter(!is.na(value))

x <- mat2 %>% select(value,ebw,ecr,stm,seo,pae,pau)
y <- x %>% select(ebw,ecr,stm,seo,pae,pau) %>% t()
RT_conn <- compare.multi.evol.rates(A=y,gp=x$value,phy=tree_nw,iter=999)




####
RT <- RT_conn
df <- data.frame(rates=RT$sigma.d.gp,
           groups=RT$groups)
ggplot(df,aes(x=as.numeric(df$groups),y=df$rates))+
  geom_point()+
  theme_minimal()+
  ggtitle(label = RT$call)+
  theme(axis.text.x = element_text(angle = 90)) 


summary(RT_conn)
plot(RT_conn)
hist(connect$value)



# Rates from tSNE clusters ~ distance
rate_tSNE_df <- data.frame(rates=RT_tSNE$sigma.d.gp,
                           clusters=as.numeric(RT_tSNE$groups))

mat_A <- full_join(mat2,clust,by="drugdrug")
mat_B <- full_join(mat_A,rate_tSNE_df,by="clusters")

rat_mat <- mat_B %>% select(Min_distance,rates,clusters) %>% 
  filter(!is.na(Min_distance))

ggplot(rat_mat,aes(x=as.factor(Min_distance),y=rates))+
  geom_point()+
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) 

rat_mat %>% group_by(Min_distance) %>% 
  summarise(mean_rate=mean(rates)) %>% 
  ggplot(aes(x=Min_distance,y=mean_rate))+
  geom_point()+
  theme_minimal()


#First, calculate the phylogenetic mean for each variable.
Sigma1 <- vcv(tree_nw) #This is the variance-covariance matrix for the tree.

gls.Ymean<-function(Y,Sigma){
  n<-length(Y)
  #Standardize sigma to reduce rounding errors
  #      tr<-sum(diag(Sigma))
  #      Sigma<-n*Sigma/tr
  #Input
  invSigma<-solve(Sigma)
  #pgls mean and variance
  X1<-rep(1,n)
  q<-2          # correct if multivariate!!
  C1<-solve(t(X1)%*%invSigma%*%X1)
  Y_PGLSmean<-C1%*%t(X1)%*%invSigma%*%Y
  Y_PGLSdeviations = Y - c(Y_PGLSmean)
  Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
  SE_Y_mean = sqrt(Y_PGLSvariance/n)
  #Save model
  results<-cbind(Y_PGLSmean,SE_Y_mean,Y_PGLSvariance)
  colnames(results)<-c("Ymean","YSE","Y_PGLSvariance")
  return(results)
}

list1 <- list()
for (i in 1:dim(dat)[2]){
  list1[[i]] <- gls.Ymean(dat[,i],Sigma=Sigma1)
}
df1 <- as.data.frame(do.call(rbind, list1))
df1$cluster <- mat2$clusters


#Calculate standard arithmetic mean across variables that belong to different clusters.
df2 <- df1 %>% group_by(cluster) %>% 
  summarise(Ymean.mean=mean(Ymean),
            YSE.mean=mean(YSE),
            Y_PGLSvariance.mean=mean(Y_PGLSvariance))

#Paired sample t tests among the means of each cluster
df1.1 <- df1 %>% filter(cluster==1) %>% select(Ymean)
df1.2 <- df1 %>% filter(cluster==2) %>% select(Ymean)

# This works for every two traits. But I have 2655
phyl.pairedttest(tree=tree_nw,x1=dat[,1],x2=dat[,2])
# We want to compare two vectors that are paired, x1 and x2. 
# Where the means in the test are the means in the clusters.
# This may be a problem, what is the tree in here? 
# Should I use a simple t-test comparing phylogenetic means across clusters?
# In this case the t-test would not be paired.

#Account for multiple testing
p.adjust(p, method = BH, n = length(p))

