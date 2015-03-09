#Code to plot diversity curves (sample-standardized and raw)
#Generates figure 1 (sample-standardized diversity), S1 (range-through comparison), and S3 (diversity by publication year)
library(RCurl)

#read subsampled diversity in all lithologies (SQS at quorum of 0.55)
all_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_div_all.csv",ssl.verifypeer = FALSE)
all_div<-read.csv(text=all_div_url)

#read subsampled diversity excluding amber (SQS at quorum of 0.55)
notamber_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_div_notamber.csv",ssl.verifypeer = FALSE)
notamber_div<-read.csv(text=notamber_div_url)

#read classical rarefaction diversity data
cr_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/CR_div.csv",ssl.verifypeer = FALSE)
cr_div<-read.csv(text=cr_div_url)

#read range-through diversity (Pull of the Recent) in all lithologies
RT_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_RTdiv_all.csv",ssl.verifypeer = FALSE)
RT_div<-read.csv(text=RT_div_url)

#read raw SIB and range-through diversity
raw_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_rawSIB_RT_div.csv",ssl.verifypeer = FALSE)
raw_div<-read.csv(text=raw_div_url)

#read Labandeira diversity data
Lab_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_Ldiv_all.csv",ssl.verifypeer = FALSE)
Labandeira_div<-read.csv(text=Lab_div_url)

#read Nicholson diversity data
Nich_div_url<-getURL("https://raw.githubusercontent.com/mclapham/insect_div/master/data/fam_Ndiv_all.csv",ssl.verifypeer = FALSE)
Nicholson_div<-read.csv(text=Nich_div_url)

#reads names and ages of time intervals
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1&limit=all")
periods<-subset(time_int,time_int$level==3)
epochs<-subset(time_int,time_int$level==4)
stages<-subset(time_int,time_int$level==5)
short_int<-rbind(epochs[1:7,],stages[23:66,])

#
#
#FIG 1: Sampling-Standardized diversity

pdf("subsamp_rich.pdf",width=8)

#plots subsampled family richness of all insects
plot(cr_div$Midpoint..Ma.,cr_div$Mean.sampled.diversity,type="n",xlim=rev(range(all_div$Midpoint..Ma.)),ylim=c(-6,max(na.omit(all_div$Mean.sampled.diversity))),xlab="Age (Ma)",ylab="Subsampled family richness")

#adds classical rarefaction data
lines(cr_div$Midpoint..Ma.,cr_div$Mean.sampled.diversity,col="gray")

#adds SQS all data
points(all_div$Midpoint..Ma.,all_div$Mean.sampled.diversity,pch=16,col="firebrick4")
lines(all_div$Midpoint..Ma.,all_div$Mean.sampled.diversity,col="firebrick4")

#adds SQS data when amber excluded
points(notamber_div$Midpoint..Ma.,notamber_div$Mean.sampled.diversity,pch=16,col="steelblue3")
lines(notamber_div$Midpoint..Ma.,notamber_div$Mean.sampled.diversity,col="steelblue3")

legend(320,125,c("All insects","Excluding amber","Classical rarefaction"),pch=c(16,16,NA),lwd=1.25,col=c("firebrick4","steelblue3","gray"),bty="n",cex=1.1)

#adds timescale with periods and smaller intervals (periods and Cenozoic epochs)
rect(periods$early_age,rep(-5,nrow(periods)),periods$late_age,rep(0,nrow(periods)),col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),-2.5,periods$abbrev)

rect(short_int$early_age,rep(-10,nrow(short_int)),short_int$late_age,rep(-5,nrow(short_int)),col=paste(short_int$color))

dev.off()
#
#
#FIG S1: Range-through diversity comparison

pdf("rt_richness_comp.pdf",width=8)

#plots range-through family richness of Nicholson data
plot(Nicholson_div$age_ma,Nicholson_div$diversity,type="o",pch=16,col="steelblue3",xlim=rev(range(RT_div$Midpoint..Ma.)),ylim=c(-20,max(Nicholson_div$diversity)),xlab="Age (Ma)",ylab="Range-through family richness")

#adds Labandeira data
points(Labandeira_div$age_ma,Labandeira_div$diversity,pch=16,col="tomato3")
lines(Labandeira_div$age_ma,Labandeira_div$diversity,col="tomato3")

#adds PBDB data
points(RT_div$Midpoint..Ma.,RT_div$Range.through.families,pch=16,col="goldenrod2")
lines(RT_div$Midpoint..Ma.,RT_div$Range.through.families,col="goldenrod2")

legend(300,700,c("PBDB","Nicholson","Labandeira"),pch=16,lwd=1.25,col=c("goldenrod2","steelblue3","tomato3"),bty="n",cex=1.25)

#adds timescale with periods and smaller intervals (periods and Cenozoic epochs)
rect(periods$early_age,rep(-20,nrow(periods)),periods$late_age,rep(0,nrow(periods)),col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),-10,periods$abbrev)

rect(short_int$early_age,rep(-40,nrow(short_int)),short_int$late_age,rep(-20,nrow(short_int)),col=paste(short_int$color))

dev.off()
#
#
#FIG S3: Diversity by publication year

#READ AND FILTER OCCURRENCES
insect_occs<-read.csv("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Insecta&show=ident,phylo,lith&limit=99999")
insect_occs<-subset(insect_occs,insect_occs$matched_rank<=5)
insect_occs<-subset(insect_occs,insect_occs$genus_reso=="" | insect_occs$genus_reso=="n. gen.")

#read references
refs<-read.csv("http://paleobiodb.org/data1.1/occs/refs.txt?base_name=Insecta&limit=all")

#file to convert intervals to 10-Myr-bins
time_url<-getURL("https://raw.githubusercontent.com/mclapham/PBDB-R-scripts/master/time_convers.csv",ssl.verifypeer = FALSE)
time_conv<-read.csv(text=time_url)

#adds myr-bin corresponding to early and late age of collection
insect_occs$myrbin_max<-time_conv$myr_bin_name[match(insect_occs$early_interval,time_conv$interval_name)]
insect_occs$myrbin_min<-time_conv$myr_bin_name[match(insect_occs$late_interval,time_conv$interval_name)]

#filters occurrences to only collections belonging to a single 10-myr-bin
filtered_occs<-subset(insect_occs,insect_occs$myrbin_max==insect_occs$myrbin_min)

#age in Ma of bin midpoints
bin_midpts<-c(327.1,312.5,302.9,60.6, 48.1, 37.1, 28.5,17.3,5.8, 140.9, 130.9, 118.7,105.8, 96.5, 88.5, 77,68, 195.6, 186.3, 177.3, 168.1, 157.8, 148.2, 294.5, 281.3, 265.8, 255.7, 249.8, 241.1, 228.2, 211.6)
names(bin_midpts)<-c("Carboniferous 3","Carboniferous 4","Carboniferous 5","Cenozoic 1","Cenozoic 2","Cenozoic 3","Cenozoic 4","Cenozoic 5","Cenozoic 6","Cretaceous 1","Cretaceous 2","Cretaceous 3","Cretaceous 4","Cretaceous 5","Cretaceous 6","Cretaceous 7","Cretaceous 8","Jurassic 1","Jurassic 2","Jurassic 3","Jurassic 4","Jurassic 5","Jurassic 6","Permian 1","Permian 2","Permian 3","Permian 4","Triassic 1","Triassic 2","Triassic 3","Triassic 4")

#adds age in Ma to each occurrence
filtered_occs$age_myr<-bin_midpts[match(filtered_occs$myrbin_max,names(bin_midpts))]

#removes occurrences without families
final_occs<-subset(filtered_occs,filtered_occs$family!="")
final_occs$family<-factor(final_occs$family)

#adds publication year to each occurrence
final_occs$pubyr<-refs$pubyr[match(final_occs$reference_no,refs$reference_no)]

#reads list of geological time periods
time_int<-read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=1")
periods<-subset(time_int,time_int$level==3)
epochs<-subset(time_int,time_int$level==4)
stages<-subset(time_int,time_int$level==5)
short_int<-rbind(epochs[1:7,],stages[23:66,])

#FUNCTIONS
#rarefaction function
rarefy<-function(occs,quota) {
  if(length(occs)>=quota) {
    length(unique(sample(occs,quota,replace=F)))
  } else return(NA)
}

#extracts occurrences published post-1900 and post-1950
post1950<-subset(final_occs,final_occs$pubyr>1950)
post1900<-subset(final_occs,final_occs$pubyr>1900)

#rarefied diversity for all occurrences
rarefied_all<-sapply(split(final_occs,final_occs$age_myr),function(x) replicate(100,rarefy(x$family,150)))
div_all<-apply(rarefied_all,2,mean)
div_sd_all<-apply(rarefied_all,2,sd)

#rarefied diversity for post-1950 occurrences
rarefied_1950<-sapply(split(post1950,post1950$age_myr),function(x) replicate(100,rarefy(x$family,150)))
div_1950<-apply(rarefied_1950,2,mean)

#proportion of occurrences published post-1950
prop_1950<-table(post1950$age_myr)/table(final_occs$age_myr)
prop_1950<-prop_1950[match(names(prop_1950),names(na.omit(div_all)))]

#rarefied diversity for post-1900 occurrences
rarefied_1900<-sapply(split(post1900,post1900$age_myr),function(x) replicate(100,rarefy(x$family,150)))
div_1900<-apply(rarefied_1900,2,mean)

#proportion of occurrences published post-1900
prop_1900<-table(post1900$age_myr)/table(final_occs$age_myr)
prop_1900<-prop_1900[match(names(prop_1900),names(na.omit(div_all)))]

#generates plot
pdf("div_pubyr.pdf",width=7)
layout(matrix(c(1,2),ncol=1),heights=c(1,3))

par(mgp=c(2,0.75,0))
par(mar=c(0,4,1,2))

plot(0,0,type="n",xaxt="n",xlab="",yaxt="n",ylab="Proportion occs",xlim=rev(range(as.numeric(names(div_all)))),ylim=c(0,1))
segments(as.numeric(names(div_1900))+1.5,0,as.numeric(names(div_1900))+1.5,prop_1900,lwd=2,col="firebrick4")
segments(as.numeric(names(div_1950))-1.5,0,as.numeric(names(div_1950))-1.5,prop_1950,lwd=2,col="steelblue3")
axis(2,at=c(0.25,0.5,0.75))
text(330,0.95,"A")

par(mar=c(4,4,0,2))

plot(names(div_all),div_all,type="o",pch=16,col="gray",xlim=rev(range(as.numeric(names(div_all)))),ylim=c(28,95),xlab="Age (Ma)",ylab="Rarefied family richness")
segments(as.numeric(names(div_all)),div_all-div_sd_all,as.numeric(names(div_all)),div_all+div_sd_all,col="gray")

points(names(div_1950),div_1950,pch=16,col="steelblue3")
lines(names(div_1950),div_1950,col="steelblue3")

points(names(div_1900),div_1900,pch=16,col="firebrick4")
lines(names(div_1900),div_1900,col="firebrick4")

legend(330,93,c("All","Post-1900 only","Post-1950 only"),lwd=1,pch=16,col=c("gray","firebrick4","steelblue3"),bty="n")
text(330,95,"B")

rect(periods$early_age,rep(29,nrow(periods)),periods$late_age,rep(32,nrow(periods)),col=paste(periods$color))
text(rowMeans(cbind(periods$early_age,periods$late_age)),30.5,periods$abbrev)
rect(short_int$early_age,rep(26,nrow(short_int)),short_int$late_age,rep(29,nrow(short_int)),col=paste(short_int$color))

dev.off()
