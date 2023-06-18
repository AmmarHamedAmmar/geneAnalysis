

library(gdata)      # for reading .xls files
library(haplo.stats)
library(adegenet)
library(stringi)
library(genetics)
library(openxlsx)
library(gdata)
library(genoPlotR)
library(readxl)


file_path <- "D:/grad/Haploview/Haploview-input_N_Ansari-Pour.xls"
haplo_Data  <- read_excel("D:/grad/Haploview/Haploview-input_N_Ansari-Pour.xls")
typeof(haplo_Data$`Sample ID`[1])

fromatted_haplo <- haplo_Data[ , 3:ncol(haplo_Data)]


geno1 <- genotype(fromatted_haplo$`1114` , sep = "")
geno2 <- genotype(fromatted_haplo$`1198` , sep = "")
geno3 <- genotype(fromatted_haplo$`1906` , sep = "")
geno4 <- genotype(fromatted_haplo$`2102` , sep = "")



haploGenoData <- makeGenotypes(data.frame(geno1 , geno2 , geno3 , geno4))

newldt <- LD( haploGenoData )
LDtable(newldt)
LDplot(newldt, distance=c(1114, 1198, 1906 , 2102))

# display the results

randomString <- stri_rand_strings(1 , 10 , pattern = "[A-Za-z0-9]")
string1 <- randomString
string2 <- "_LD_plot.pdf"
fileName <- paste(string1, string2, sep ="")


# start visualizing  and store this result in pdf 
pdf(fileName)

LDtable(newldt)                            # graphical color-coded table
LDplot(newldt, distance=c(1114, 1198, 1906 , 2102))  # LD plot vs distance



# more markers makes prettier plots!
data <- list()
nobs <- 1000
ngene <- 20
s <- seq(0,1,length=ngene)
a1 <- a2 <- matrix("", nrow=nobs, ncol=ngene)
for(i in 1:length(s) )
{
  
  rallele <- function(p) sample( c("T","G"), 1, p=c(p, 1-p))
  
  if(i==1)
  {
    a1[,i] <- sample( c("T","G"), 1000, p=c(0.5,0.5), replace=TRUE)
    a2[,i] <- sample( c("T","G"), 1000, p=c(0.5,0.5), replace=TRUE)
  }
  else
  {
    p1 <- pmax( pmin( 0.25 + s[i] * as.numeric(a1[,i-1]=="T"),1 ), 0 )
    p2 <- pmax( pmin( 0.25 + s[i] * as.numeric(a2[,i-1]=="T"),1 ), 0 )
    a1[,i] <- sapply(p1, rallele )
    a2[,i] <- sapply(p2, rallele )
  }
  
  data[[paste("C",i,sep="")]] <- genotype(a1[,i],a2[,i])
}
data <- data.frame(data)
data <- makeGenotypes(data)

ldt <- LD(data)
plot(ldt, digits=2, marker=19) # d



parent_Title = "Pairwise LD -----------"


title1 = "LD Test : "
N =  paste0("N is : " , unlist(ldt$n))
P_value = paste0("P-value is : " , unlist(ldt$`P-value`))
X_2 = paste0("X^2 is : " , unlist(ldt$`X^2`))


title2 = "Estimates : "
D1 = paste0("D` is : " , unlist(ldt$`D'`))
D2 =paste0("D is : " , unlist(ldt$D))

plot(NA, xlim=c(0,9), ylim=c(0,9), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')


text(1,8 , parent_Title , pos = 4 )


text(1,7,title1 , pos = 4 )
text(1,6,N, pos=4)
text(1,5,P_value, pos=4)
text(1,4,X_2, pos=4)

text(1,3,title2 , pos = 4 ) 
text(1,2,D1, pos=4)
text(1,1,D2, pos=4)

points(rep(1,8),1:8, pch=8)


dev.off()











