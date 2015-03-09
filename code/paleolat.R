#Code to analyze paleolatitudinal distribution of insect fossil occurrences
#plots figure 3 (paleolatitude histograms) and figure S2 (paleolatitude of occurrences through time)
library(RCurl)

#DATA ACQUISITION
#reads insect occurrences
insect_occs<-read.csv("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Insecta&show=phylo,loc,paleoloc,lith,geo&limit=99999")

#reads file to map time intervals to 10-Ma bins
time_url<-getURL("https://raw.githubusercontent.com/mclapham/PBDB-R-scripts/master/time_convers.csv",ssl.verifypeer = FALSE)
time_conv<-read.csv(text=time_url)

#DATA PROCESSING
#converts time intervals to 10-Ma bins
insect_occs$myrbin_max<-time_conv$myr_bin_name[match(insect_occs$early_interval,time_conv$interval_name)]
insect_occs$myrbin_min<-time_conv$myr_bin_name[match(insect_occs$late_interval,time_conv$interval_name)]

#filters to exclude occurrences not resolved to a single bin
filtered_occs<-subset(insect_occs,insect_occs$myrbin_max==insect_occs$myrbin_min)

#10-Ma bin midpoint ages
bin_midpts<-c(327.1,312.5,302.9,60.6, 48.1, 37.1, 28.5,17.3,5.8, 140.9, 130.9, 118.7,105.8, 96.5, 88.5, 77,68, 195.6, 186.3, 177.3, 168.1, 157.8, 148.2, 294.5, 281.3, 265.8, 255.7, 249.8, 241.1, 228.2, 211.6)
names(bin_midpts)<-c("Carboniferous 3","Carboniferous 4","Carboniferous 5","Cenozoic 1","Cenozoic 2","Cenozoic 3","Cenozoic 4","Cenozoic 5","Cenozoic 6","Cretaceous 1","Cretaceous 2","Cretaceous 3","Cretaceous 4","Cretaceous 5","Cretaceous 6","Cretaceous 7","Cretaceous 8","Jurassic 1","Jurassic 2","Jurassic 3","Jurassic 4","Jurassic 5","Jurassic 6","Permian 1","Permian 2","Permian 3","Permian 4","Triassic 1","Triassic 2","Triassic 3","Triassic 4")

#adds numeric age to occurrences
filtered_occs$age_myr<-bin_midpts[match(filtered_occs$myrbin_max,names(bin_midpts))]

#extracts compression/impression occurrences
no_amber<-subset(filtered_occs,filtered_occs$lithology1!="amber")

#subsets occurrences in northern hemisphere regions
fmr_soviet_lat<-subset(no_amber$paleolat,no_amber$cc %in% c("RU","KZ","KG","TJ","UZ","UA","BY","AZ","AM","GE","LT","LV","EE","MD"))
n_am_lat<-subset(no_amber$paleolat,no_amber$cc %in% c("CA","US","MX"))
europe_lat<-subset(no_amber$paleolat,no_amber$cc %in% c("UK","FR","DE","NL","BE","ES","PT","DK","CH","LU","AT","LI","CZ","IT","PL","HU","RO","HR","SK","SI"))

#finds min and max paleolatitude of entire data
max_lat<-max(c(fmr_soviet_lat,n_am_lat,europe_lat),na.rm=T)
min_lat<-min(c(fmr_soviet_lat,n_am_lat,europe_lat),na.rm=T)

#plots figure 3
pdf("paleolatitude.pdf",width=5,height=8)
par(mfrow=c(3,1))
par(mar=c(4.5,4,2,0))
par(mgp=c(2,0.75,0))
hist(n_am_lat,xlim=c(min_lat,max_lat),breaks=20,col="gray",xlab="Paleolatitude",ylab="Number of occurrences",main="North America")
hist(europe_lat,xlim=c(min_lat,max_lat),breaks=20,col="gray",xlab="Paleolatitude",ylab="Number of occurrences",main="Europe")
hist(fmr_soviet_lat,xlim=c(min_lat,max_lat),breaks=20,col="gray",xlab="Paleolatitude",ylab="Number of occurrences",main="Former Soviet Union")
dev.off()

#plots figure S2
pdf("median_paleolat.pdf",width=8,height=8)
par(mgp=c(2,0.75,0))
plot(filtered_occs$age_myr,filtered_occs$paleolat,pch=16,col=rgb(0.75,0.75,0.75,0.1),xlim=rev(range(filtered_occs$age_myr,na.rm=T)),xlab="Age (Ma)",ylab="Paleolatitude")
median_lat <- sapply(split(filtered_occs$paleolat,filtered_occs$age_myr),function(x) median(x,na.rm=T))
points(names(median_lat),median_lat,pch=15)
dev.off()
