library(RCurl)

#READ AND FILTER OCCURRENCES
insect_occs<-read.csv("http://paleobiodb.org/data1.1/occs/list.txt?base_name=Insecta&show=ident,phylo,lith&limit=99999")
insect_occs<-subset(insect_occs,insect_occs$matched_rank=="species")
insect_occs<-subset(insect_occs,insect_occs$genus_reso=="" | insect_occs$genus_reso=="n. gen.")

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
filtered_occs<-subset(filtered_occs,filtered_occs$family!="")

mean_sp<-numeric(length(bin_midpts))
occ_ct<-numeric(length(bin_midpts))

pdf("species_per_fam.pdf",width=8)
plot(1,1,type="n",log="y",xlab="Age (Ma)",ylab="Species per family",xlim=rev(range(bin_midpts)),ylim=c(1,250))

for (i in 1:length(bin_midpts)) {
  temp_list<-subset(filtered_occs,filtered_occs$age_myr==sort(bin_midpts)[i])
  temp_list$family<-factor(temp_list$family)
  sp_per_fam<-sapply(split(temp_list$matched_name,factor(temp_list$family)),function(y) length(unique(y)))
  points(rep(sort(bin_midpts)[i],length(sp_per_fam)),sp_per_fam,pch=21,bg=rgb(0,0,0,0.05))
  points(sort(bin_midpts)[i],mean(sp_per_fam),pch=16,col="red",cex=1.5)
  mean_sp[i]<-mean(sp_per_fam)
  occ_ct[i]<-nrow(temp_list)
}
dev.off()

#plot fig S5 (species per family data)
pdf("sp_fam_correlation.pdf",width=6)
par(mgp=c(2,0.75,0))
plot(occ_ct,mean_sp,log="x",xlim=c(100,6500),xlab="Occurrences",ylab="Mean species count per family",pch=16)
dev.off()
