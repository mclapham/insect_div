#Code to analyze relationship between tectonic setting and SQS subsampled diversity
#Plots figure 2
library(RCurl)

#DATA ACQUISITION
#reads insect occurrences
insect_occs<-read.csv("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Insecta&show=phylo,loc,paleoloc,lith,geo&limit=99999")

#reads file to map time intervals to 10-Ma bins
time_url<-getURL("https://raw.githubusercontent.com/mclapham/PBDB-R-scripts/master/time_convers.csv",ssl.verifypeer = FALSE)
time_conv<-read.csv(text=time_url)

#reads file with additional tectonic settings
tect_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/tectonic_settings.csv",ssl.verifypeer = FALSE)
tectonic<-read.csv(text=tect_url)

#reads sampling-standardized diversity results from compression/impression fossils
no_amber_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_div_notamber.csv",ssl.verifypeer = FALSE)
no_amber_div<-read.csv(text=no_amber_url)


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

#adds tectonic settings to PBDB collections missing that data (because I lack editing permission)
filtered_occs$tectonic_setting[which(filtered_occs$collection_no %in% tectonic$collection_no)]<-tectonic$tectonic_setting[na.omit(match(filtered_occs$collection_no,tectonic$collection_no))]

#extracts compression/impression occurrences
no_amber<-subset(filtered_occs,filtered_occs$lithology1!="amber")


#DATA ANALYSIS
#counts occurrences in "high-subsidence" basin types
hsb_count<-sapply(split(no_amber$tectonic_setting,no_amber$age_myr),function(x) length(subset(x,x %in% c("rift","volcanic basin","pull-apart basin"))))

#counts occurrences in "low-subsidence" basin types
lsb_count<-sapply(split(no_amber$tectonic_setting,no_amber$age_myr),function(x) length(x)-length(subset(x,x %in% c("rift","volcanic basin","pull-apart basin"))))

#first differences (change from one time interval to succeeding interval)
hsb_change<-diff(rev(hsb_count)) #change in high-subsidence basin occurrence counts
lsb_change<-diff(rev(lsb_count)) #change in low-subsidence basin occurrence counts
diversity_ch<-diff(rev(no_amber_div$Mean.sampled.diversity)) #change in SQS diversity

#linear regression
summary(lm(diversity_ch ~ hsb_change + lsb_change))

#plot change in diversity vs. change in occurrences in high- and low-subsidence basins
pdf("diversity_tectonics.pdf",width=5,height=9)
par(mfrow=c(2,1))
par(mar=c(4,4,1,2))
par(mgp=c(2,0.75,0))
plot(hsb_change,diversity_ch,pch=16,cex=1.25,cex.axis=0.75,xlab=expression(paste(Delta," occurrences in high-subsidence basins")),ylab=expression(paste(Delta," subsampled family diversity")))
plot(lsb_change,diversity_ch,pch=16,cex=1.25,cex.axis=0.75,xlab=expression(paste(Delta," occurrences in low-subsidence basins")),ylab=expression(paste(Delta," subsampled family diversity")))
dev.off()
