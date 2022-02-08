
source("Urban-Rural Parameter Inference.R")

#Assign colors
# Function to assign transparent colors in base R
makeTransparent <- function(someColor, alpha = 100){
  if(alpha < 1){alpha <- alpha*100} else{alpha <- alpha} # always forget if alpha is a decimal or whole number
  ccdfun <- function(x){rgb(red = x[1], green = x[2], blue = x[3], alpha = alpha, maxColorValue = 255)}
  apply(col2rgb(someColor), 2, ccdfun)
}

colors <- c("#ADB6B6E5", "#AD002AE5", "#925E9FE5")
col_total <- colors[1]
col_urban <- colors[2]
col_rural <- colors[3]
colt_total <- makeTransparent(col_total, alpha = 60)
colt_urban <- makeTransparent(col_urban, alpha = 60)
colt_rural <- makeTransparent(col_rural, alpha = 60)

#####################################################################
## Absolute probability of cluster of size 30
#####################################################################
# Function to calculate probability of a cluster of size Y
clustprob <- function(y, R, k){
  l <- lgamma(k*y + (y-1)) - (lgamma(k*y) + lgamma(y+1)) + (y-1) * log(R/k) - (k*y + (y-1)) * log(1 + R/k)
  return(exp(l))
}

#Find max cluster size in the data
maxobservedclust <- max(total[,1])

# Calculate probabilities for each cluster size from 1:maxclustersize
total.probclust <- urban.probclust <- rural.probclust <- matrix(NA, nrow = maxobservedclust, ncol = 3)
for (j in 1:maxobservedclust){
  urban.probclust[j,1] <- clustprob(j, urban.ests[1,1], urban.ests[2,1])
  urban.probclust[j,2] <- clustprob(j, urban.ests[1,2], urban.ests[2,2])
  urban.probclust[j,3] <- clustprob(j, urban.ests[1,3], urban.ests[2,3])
  rural.probclust[j,1] <- clustprob(j, rural.ests[1,1], rural.ests[2,1])
  rural.probclust[j,2] <- clustprob(j, rural.ests[1,2], rural.ests[2,2])
  rural.probclust[j,3] <- clustprob(j, rural.ests[1,3], rural.ests[2,3])
}

maxrange <- seq(1, maxobservedclust, 1)
xx <- seq(0, maxobservedclust, by = 0.01)

# Probability of a cluster of size Y
plot(log(range(1, maxobservedclust)), range(-11,0), type = 'n', xlab = '', ylab = '', axes = F)
  points(log(maxrange), log(urban.probclust[,1]), col = col_urban, pch = 15, cex = 1)
  points(log(maxrange), log(rural.probclust[,1]), col = col_rural, pch = 17, cex = 1)
  axis(side = 1, at = log(c(1, 5, seq(10,maxobservedclust+10,b=10))), c(1, 5, seq(10, maxobservedclust+10,b=10)))
  axis(side = 2, at = log(c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1)), c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1))
  mtext(side = 1, 'Total Outbreak Size', padj=4)
  mtext(side = 2, 'Probability', padj=-4)
  legend('topright', c("Urban","Rural"), 
         pch=c(15, 17),cex=1, col = c(col_urban, col_rural), bty="n")

# underlying gamma distributed nu values
plot(log(1+range(0,maxobservedclust)), range(-12,0), type='n', axes=F, xlab='', ylab='')
# urban model
  lines(log(xx),log(dgamma(xx, scale = urban.ests[1,1]/urban.ests[2,1], shape = urban.ests[2,1])), lwd = 2, lty=1, col=col_urban)
  polygon(c(log(xx), rev(log(xx))),c(log(dgamma(xx, scale = urban.ests[1,2]/urban.ests[2,2], shape = urban.ests[2,2])),rev(log(dgamma(xx, scale = urban.ests[1,3]/urban.ests[2,3], shape = urban.ests[2,3])))),col = colt_urban, border=NA)
# rural model
  lines(log(xx),log(dgamma(xx, scale = rural.ests[1,1]/rural.ests[2,1], shape = rural.ests[2,1])), lwd = 2,lty=1, col=col_rural)
  polygon(c(log(xx), rev(log(xx))), c(log(dgamma(xx, scale = rural.ests[1,2]/rural.ests[2,1], shape = rural.ests[2,1])),rev(log(dgamma(xx, scale = rural.ests[1,3]/rural.ests[2,1], shape = rural.ests[2,1])))),col = colt_rural, border=NA)
  axis(side=1, at=log(c(1, 5, seq(10,maxobservedclust+10,b=10))), c(1, 5, seq(10, maxobservedclust+10,b=10)))
  axis(side=2, at=log(c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1)), c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1))
  mtext(side=1, expression('Individual Reproductive Number,' ~ nu), padj = 4)
  mtext(side=2, 'Density', padj = -4)
  legend('topright', c("Urban","Rural"), lty=1, lwd=2, cex=1, col=c(col_urban, col_rural), bty="n")
  legend('topright', c(rep(NA,2),"Uncertainty\nInterval"), fill=c(rep(NA,2),"grey"), bty="n", cex=1, border=c(rep(NA,2),"black"))

#############################
  
L.set <- seq(2,50,1)
p.urban <- p.rural <- p.ur.rr <-  c()
p.urban.lo <- p.rural.lo <- p.ur.rr.lo <- c()
p.urban.hi <- p.rural.hi <- p.ur.rr.hi <- c()
  
for (i in 1:length(L.set)){
    p.urban[i] <-  (1 - sum(clustprob(1:(L.set[i]-1), urban.ests[1,1], urban.ests[2,1]))) # MLE
    p.rural[i] <-  (1 - sum(clustprob(1:(L.set[i]-1), rural.ests[1,1], rural.ests[2,1]))) # MLE
    p.ur.rr[i] <- p.rural[i] / p.urban[i] # Relative
    
    p.urban.lo[i] <-  (1 - sum(clustprob(1:(L.set[i]-1), urban.ests[1,2], urban.ests[2,2]))) # Lowest R and k
    p.rural.lo[i] <-  (1 - sum(clustprob(1:(L.set[i]-1), rural.ests[1,2], rural.ests[2,2]))) # Lowest R and k
    p.ur.rr.lo[i] <- p.rural.lo[i] / p.urban.lo[i] # Relative
    
    p.urban.hi[i] <-  (1 - sum(clustprob(1:(L.set[i]-1),urban.ests[1,3],urban.ests[2,3]))) # Highest R and k
    p.rural.hi[i] <-  (1 - sum(clustprob(1:(L.set[i]-1),rural.ests[1,3],rural.ests[2,3]))) # highest R and k
    p.ur.rr.hi[i] <- p.rural.hi[i] / p.urban.hi[i] # Relative
  }
  
# Relative probability
plot((range(1, max(L.set))), range(0,50), type='n', xlab='', ylab='', axes=F)
  lines((L.set), (p.ur.rr), col = col_urban, lwd=3, lty=1)
  polygon(c(L.set, rev(L.set)), c(p.ur.rr.lo, rev(p.ur.rr.hi)), col = colt_urban, border = NA)
  lines((L.set), (p.ur.rr.lo), col = col_urban, lwd = 1, lty = 1)
  lines((L.set), (p.ur.rr.hi), col = col_urban, lwd = 1, lty = 1)
  axis(side = 1); axis(side = 2)
  mtext(side = 1, 'Minimum Total Outbreak Size', padj=4)
  mtext(side = 2, 'Relative Probability ', padj=-4)
  legend('topleft', "MLE estimate", lty = 1, lwd = 3, col = col_urban, bty='n')
  legend('topleft', c(NA,"    Uncertainty Interval"), fill=c(NA,colt_urban), bty="n", cex=1, border=c(NA,col_urban))

#####################################################################
## Absolute probability of cliuster of size 30
#####################################################################
L <- 30
Rs <- seq(0.02,1.0,0.005)
ks <-  10^seq(-2.3,0,0.005)
P30 <- matrix(NA, ncol = length(ks), nrow = length(Rs))

# This calculates the probability of observing a cluster of size 30 for all values across Rs and ks
for(i in 1:length(Rs)) {
  for(j in 1:length(ks)) {
    P30[i,j] <- 1 - sum(clustprob(1:(L-1),Rs[i],ks[j]))
  }
}

# Set up figure 
purb <- round(1 - sum(clustprob(1:(L-1),urban.ests[1,1],urban.ests[2,1])),4)
prur <- round(1 - sum(clustprob(1:(L-1),rural.ests[1,1],rural.ests[2,1])),4)
lbreaks <- c(-1e-10,1e-5,1e-3,5e-3,1e-2,2e-2,5e-2,1e-1,1.5e-1,1)

image(Rs, ks, P30, xlab = bquote("Reproduction Number," ~ italic("R")), ylab = bquote("Dispersion parameter," ~ italic("k")),
      log="y",col = rev(heat.colors(length(lbreaks)-1)), breaks = lbreaks)
  contour(Rs, ks, P30, levels = lbreaks, labcex = 0.8, add = TRUE)
  points(urban.ests[1,1],urban.ests[2,1], pch = 21, col = "black", bg = col_urban, cex = 2)
  points(rural.ests[1,1],rural.ests[2,1], pch = 21, col = "black", bg = col_rural, cex = 2)
  text(urban.ests[1,1], urban.ests[2,1], purb)
  text(rural.ests[1,1], rural.ests[2,1], prur)
  legend('topleft', c("Urban","Rural"), 
       pch = 21, col = "black", pt.bg = c(col_urban,col_rural), bty = 'n', cex = 1.5)


#####################################################################
## Number of introductions until first large outbreak  
#####################################################################
bp <- function(gens = Inf, init.size = 1, offspring, ...){  
  Z <- list() #initiate the list
  Z[[1]] <- init.size #set the first position of the list as the number of index cases
  i <- 1 
  while(sum(Z[[i]]) > 0 && i <= gens) { 
    Z[[i+1]] <- offspring(sum(Z[[i]]), ...) 
    i <- i+1 
  } 
  return(Z)
} 
  
num_sims <- 500
num_chains <- 2000

# Identify the number of introductions until first cluster of size 15, 30, and 50
set.seed(05062020)
max_clust <- c(15, 30, 50)

first.sse.urban.15 <- first.sse.rural.15 <- c()

for (i in 1:num_sims){
  z.urban <- replicate(num_chains, bp(offspring = rnbinom, mu = urban.ests[1,1], size = urban.ests[2,1])) # Simulate surveillance system
  y.urban <- unlist(lapply(z.urban, function(x) sum(unlist(x)))) # generate cluster data
  sse.u <- cbind(seq_along(y.urban), ifelse(y.urban > max_clust[1], 1, 0)) # indicate if the cluster size is greater than the specified max_clust[1]er (here, 15)
  z.rural <- replicate(num_chains, bp(offspring = rnbinom, mu = rural.ests[1,1], size = rural.ests[2,1]))
  y.rural <- unlist(lapply(z.rural, function(x) sum(unlist(x))))
  sse.r <- cbind(seq_along(y.rural), ifelse(y.rural > max_clust[1], 1, 0))
  first.sse.urban.15[i] <- min(sse.u[sse.u[,2] == 1, 1]) # Identify the position of the first SSE (number of introductions to first SSE)
  first.sse.rural.15[i] <- min(sse.r[sse.r[,2] == 1, 1])
}

first.sse.urban.30 <- first.sse.rural.30 <- c()
for (i in 1:num_sims){
  z.urban <- replicate(num_chains, bp(offspring = rnbinom,mu = urban.ests[1,1], size = urban.ests[2,1]))
  y.urban <- unlist(lapply(z.urban, function(x) sum(unlist(x))))
  sse.u <- cbind(seq_along(y.urban), ifelse(y.urban > max_clust[2], 1, 0))
  
  z.rural <- replicate(num_chains, bp(offspring = rnbinom, mu = rural.ests[1,1], size = rural.ests[2,1]))
  y.rural <- unlist(lapply(z.rural, function(x) sum(unlist(x))))
  sse.r <- cbind(seq_along(y.rural), ifelse(y.rural > max_clust[2], 1, 0))
  
  first.sse.urban.30[i] <- min(sse.u[sse.u[,2] == 1,1])
  first.sse.rural.30[i] <- min(sse.r[sse.r[,2] == 1,1])
}

first.sse.urban.50 <- first.sse.rural.50 <- c()
for (i in 1:num_sims){
  z.urban <- replicate(num_chains, bp(offspring = rnbinom, mu = urban.ests[1,1], size = urban.ests[2,1]))
  y.urban <- unlist(lapply(z.urban, function(x) sum(unlist(x))))
  sse.u <- cbind(seq_along(y.urban),ifelse(y.urban > max_clust[3], 1, 0))
  z.rural <- replicate(num_chains, bp(offspring = rnbinom, mu = rural.ests[1,1], size = rural.ests[2,1]))
  y.rural <- unlist(lapply(z.rural, function(x) sum(unlist(x))))
  sse.r <- cbind(seq_along(y.rural), ifelse(y.rural > max_clust[3], 1, 0))
  first.sse.urban.50[i] <- min(sse.u[sse.u[,2] == 1, 1])
  first.sse.rural.50[i] <- min(sse.r[sse.r[,2] == 1, 1])
}

# Organize the data
urban.15 <- data.frame(first.sse=first.sse.urban.15[first.sse.urban.15 < num_chains], 
                       model="Urban Model", 
                       sse_size="Outbreak Size, Y=15")
rural.15 <- data.frame(first.sse=first.sse.rural.15[first.sse.rural.15 < num_chains], 
                       model="Rural Model", 
                       see_size="Outbreak Size, Y=15")
urban.30 <- data.frame(first.sse=first.sse.urban.30[first.sse.urban.30 < num_chains], 
                       model="Urban Model", 
                       sse_size="Outbreak Size, Y=30")
rural.30 <- data.frame(first.sse=first.sse.rural.30[first.sse.rural.30 < num_chains], 
                       model="Rural Model", 
                       see_size="Outbreak Size, Y=30")
urban.50 <- data.frame(first.sse=first.sse.urban.50[first.sse.urban.50 < num_chains], 
                       model="Urban Model", 
                       sse_size="Outbreak Size, Y=50")
rural.50 <- data.frame(first.sse=first.sse.rural.50[first.sse.rural.50 < num_chains], 
                       model="Rural Model", see_size="Outbreak Size, Y=50")

nnames <- c("first.sse", "Model","sse_size")
colnames(urban.15) <- colnames(urban.30) <- colnames(urban.50) <- nnames
colnames(rural.15) <- colnames(rural.30) <- colnames(rural.50) <- nnames

violindata <- rbind(urban.15, urban.30, urban.50, rural.15, rural.30, rural.50)

library(ggplot2)
ggplot(data = violindata, aes(Model, first.sse, fill = Model)) + 
  theme_classic() + theme(legend.position = "none") +
  ylab("Number of Incident Cases Until First Large Outbreak") +
  geom_violin(alpha = 0.8) + 
  geom_boxplot(width = 0.05, fill ="white")+
  coord_flip() +
  facet_grid(~sse_size)+
  scale_fill_manual(values=c(col_rural, col_urban))


