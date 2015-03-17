#REQUIRED SETTINGS AND PACKAGES
setInternet2(use=FALSE)
library(RCurl)

#READ AND FILTER OCCURRENCES
insect_occs <- read.csv("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Insecta&show=ident,phylo,lith&limit=99999")
insect_occs <- subset(insect_occs,insect_occs$matched_rank %in% c("genus","subgenus","species"))
insect_occs <- subset(insect_occs,insect_occs$genus_reso=="" | insect_occs$genus_reso=="n. gen.")

#file to convert intervals to 10-Myr-bins
time_url <- getURL("https://raw.githubusercontent.com/mclapham/PBDB-R-scripts/master/time_convers.csv",ssl.verifypeer = FALSE)
time_conv <- read.csv(text=time_url)

#adds myr-bin corresponding to early and late age of collection
insect_occs$myrbin_max <- time_conv$myr_bin_name[match(insect_occs$early_interval,time_conv$interval_name)]
insect_occs$myrbin_min <- time_conv$myr_bin_name[match(insect_occs$late_interval,time_conv$interval_name)]

#filters occurrences to only collections belonging to a single 10-myr-bin
filtered_occs <- subset(insect_occs,insect_occs$myrbin_max==insect_occs$myrbin_min)

#age in Ma of bin midpoints
bin_midpts <- c(327.1,312.5,302.9,60.6, 48.1, 37.1, 28.5,17.3,5.8, 140.9, 130.9, 118.7,105.8, 96.5, 88.5, 77,68, 195.6, 186.3, 177.3, 168.1, 157.8, 148.2, 294.5, 281.3, 265.8, 255.7, 249.8, 241.1, 228.2, 211.6)
names(bin_midpts) <- c("Carboniferous 3","Carboniferous 4","Carboniferous 5","Cenozoic 1","Cenozoic 2","Cenozoic 3","Cenozoic 4","Cenozoic 5","Cenozoic 6","Cretaceous 1","Cretaceous 2","Cretaceous 3","Cretaceous 4","Cretaceous 5","Cretaceous 6","Cretaceous 7","Cretaceous 8","Jurassic 1","Jurassic 2","Jurassic 3","Jurassic 4","Jurassic 5","Jurassic 6","Permian 1","Permian 2","Permian 3","Permian 4","Triassic 1","Triassic 2","Triassic 3","Triassic 4")

#adds age in Ma to each occurrence
filtered_occs$age_myr <- bin_midpts[match(filtered_occs$myrbin_max,names(bin_midpts))]

final_occs <- subset(filtered_occs,filtered_occs$family!="")
final_occs$family <- factor(final_occs$family)
                          
#reads list of geological time periods
time_int <- read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1")
periods <- subset(time_int,time_int$level==3)
epochs <- subset(time_int,time_int$level==4)
stages <- subset(time_int,time_int$level==5)
short_int <- rbind(epochs[1:7,],stages[23:66,])

#similarity function
forbes <- function(x,y)  {
  a <- length(which((x * y)>0))
  b <- length(which(x>0))-a
  c <- length(which(y>0))-a
  n <- a+b+c
  return(a*n/((a + b)*(a + c)))
}

occs_matrix <- sapply(split(final_occs$family,final_occs$age_myr),table)
occs_matrix[which(occs_matrix>0)] <- 1 #sets all counts to 1 for presence

#creates variable for similarity values
bin_simil <- numeric(ncol(occs_matrix)-2)

#loops through matrix to calculate similarity between time interval and preceding
for (i in 1:(ncol(occs_matrix)-2)) {
  bin_simil[i] <- forbes(occs_matrix[,i], occs_matrix[,i+1])
}

#Ages of bins
bin_ages <- as.numeric(colnames(occs_matrix)[1:(ncol(occs_matrix)-2)])

#Figure S4
pdf("bin_similarity.pdf",width=6)

par(mgp=c(2,0.75,0))
plot(bin_ages,bin_simil,xlab="Age (Ma)",ylab="Similarity to preceding bin",xlim=rev(range(bin_ages)),ylim=c(0.43,0.9))

rect(periods$early_age,rep(0.44,nrow(periods)),periods$late_age,rep(0.46,nrow(periods)),col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),0.45,periods$abbrev)
rect(short_int$early_age,rep(0.42,nrow(short_int)),short_int$late_age,rep(0.44,nrow(short_int)),col=paste(short_int$color))
dev.off()

