#this script plot one kmer at a time, can plot intact k genome hits or splitted k genome hits

setwd("/Users/dorothytam/Documents/Bath_postdoc/genomic_analyses/170_USAgenomes_analysis/111_yearGWAS//")
mypheno<-read.table("phenotype.txt",header=T)
table(mypheno$year)
#0  1 
#44 67

mytable<-read.table("111yearGWAS_mystartendLR_k200.txt",header=T)
mymerge<-merge(mytable,mypheno,by.x="genome",by.y="ID")
colnames(mymerge)[8]<-"case_control"

hist(mytable$EndR)

unique(mytable$kmer)[1:50]
[1] "kmer1"   "kmer2"   "kmer3"   "kmer4"   "kmer5"   "kmer7"  
[7] "kmer9"   "kmer11"  "kmer12"  "kmer13"  "kmer15"  "kmer17" 
[13] "kmer18"  "kmer19"  "kmer21"  "kmer27"  "kmer28"  "kmer31" 
[19] "kmer32"  "kmer33"  "kmer34"  "kmer35"  "kmer36"  "kmer37" 
[25] "kmer39"  "kmer42"  "kmer43"  "kmer44"  "kmer45"  "kmer46" 
[31] "kmer47"  "kmer50"  "kmer51"  "kmer52"  "kmer53"  "kmer56" 
[37] "kmer81"  "kmer83"  "kmer85"  "kmer86"  "kmer87"  "kmer88" 
[43] "kmer89"  "kmer90"  "kmer91"  "kmer92"  "kmer94"  "kmer95" 
[49] "kmer98"  "kmer100"

mykmer="kmer93"  #pick the kmer to plot

mykmer_tab<-mymerge[which(mymerge$kmer==mykmer),]
myintactk<-mykmer_tab[which(abs(mykmer_tab$EndL-mykmer_tab$StartR)<10000),]
mysplitk<-mykmer_tab[which(abs(mykmer_tab$EndL-mykmer_tab$StartR)>10000),]

#pick to plot the kmer genome-hit with splitted flanks or intact k
#myflk_behave_pheno<-myintactk  
myflk_behave_pheno<-mysplitk

###preparing the plotting coordinates for ctrlL, ctrlR, caseL and caseR###

myktab_ctrl<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer & myflk_behave_pheno$case_control=="0"),]
myktab_case<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer & myflk_behave_pheno$case_control=="1"),]

#processing the left flank, make a table storing the coordinate for each arrow and the size
myktab<-myktab_ctrl[order(myktab_ctrl$StartL,decreasing=F),]
breakpt<-which(diff(myktab$StartL)>10000)

if(length(breakpt)>0){ #if at least one break point found
  mymed_ctrl_L<-matrix(0,length(breakpt)+1,3)
  colnames(mymed_ctrl_L)<-c("startmed","endmed","count")
  mybreak<-c(1,breakpt,nrow(myktab))
  for(j in 1:(length(mybreak)-1)){
    print(j)
    if(j==1){
      mySbreak<-1
      myEbreak<-mybreak[j+1]
    }else { #if the break value is not the first nor the last
      mySbreak<-mybreak[j]+1
      myEbreak<-mybreak[j+1]
    }
    mystart_med<-median(myktab$StartL[mySbreak:myEbreak])
    myend_med<-median(myktab$EndL[mySbreak:myEbreak])
    mymed_ctrl_L[j,"startmed"]<-mystart_med
    mymed_ctrl_L[j,"endmed"]<-myend_med
    mymed_ctrl_L[j,"count"]<-length(mySbreak:myEbreak)
  }
}

if(length(breakpt)==0){    #if no break point found
  mymed_ctrl_L<-matrix(0,1,3)
  colnames(mymed_ctrl_L)<-c("startmed","endmed","count")
  mySbreak<-1
  myEbreak<-nrow(myktab)
  mystart_med<-median(myktab$StartL[mySbreak:myEbreak])
  myend_med<-median(myktab$EndL[mySbreak:myEbreak])
  mymed_ctrl_L[1,"startmed"]<-mystart_med
  mymed_ctrl_L[1,"endmed"]<-myend_med
  mymed_ctrl_L[1,"count"]<-length(mySbreak:myEbreak)
}

#processing the right flank, make a table storing the coordinate for each arrow and the size
myktab<-myktab_ctrl[order(myktab_ctrl$StartR,decreasing=F),]
breakpt<-which(diff(myktab$StartR)>10000)

if(length(breakpt)>0){ #if at least one break point found
  mymed_ctrl_R<-matrix(0,length(breakpt)+1,3)
  colnames(mymed_ctrl_R)<-c("startmed","endmed","count")
  mybreak<-c(1,breakpt,nrow(myktab))
  for(j in 1:(length(mybreak)-1)){
    print(j)
    if(j==1){
      mySbreak<-1
      myEbreak<-mybreak[j+1]
    }else { #if the break value is not the first nor the last
      mySbreak<-mybreak[j]+1
      myEbreak<-mybreak[j+1]
    }
    mystart_med<-median(myktab$StartR[mySbreak:myEbreak])
    myend_med<-median(myktab$EndR[mySbreak:myEbreak])
    mymed_ctrl_R[j,"startmed"]<-mystart_med
    mymed_ctrl_R[j,"endmed"]<-myend_med
    mymed_ctrl_R[j,"count"]<-length(mySbreak:myEbreak)
  }
}

if(length(breakpt)==0){    #if no break point found
  mymed_ctrl_R<-matrix(0,1,3)
  colnames(mymed_ctrl_R)<-c("startmed","endmed","count")
  mySbreak<-1
  myEbreak<-nrow(myktab)
  mystart_med<-median(myktab$StartR[mySbreak:myEbreak])
  myend_med<-median(myktab$EndR[mySbreak:myEbreak])
  mymed_ctrl_R[1,"startmed"]<-mystart_med
  mymed_ctrl_R[1,"endmed"]<-myend_med
  mymed_ctrl_R[1,"count"]<-length(mySbreak:myEbreak)
}

#extract the table of the kmer and case
myktab_case<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer & myflk_behave_pheno$case_control=="1"),]

#processing the left flank, make a table storing the coordinate for each arrow and the size
myktab<-myktab_case[order(myktab_case$StartL,decreasing=F),]
breakpt<-which(diff(myktab$StartL)>10000)

if(length(breakpt)>0){ #if at least one break point found
  mymed_case_L<-matrix(0,length(breakpt)+1,3)
  colnames(mymed_case_L)<-c("startmed","endmed","count")
  mybreak<-c(1,breakpt,nrow(myktab))
  for(j in 1:(length(mybreak)-1)){
    print(j)
    if(j==1){
      mySbreak<-1
      myEbreak<-mybreak[j+1]
    }else { #if the break value is not the first nor the last
      mySbreak<-mybreak[j]+1
      myEbreak<-mybreak[j+1]
    }
    mystart_med<-median(myktab$StartL[mySbreak:myEbreak])
    myend_med<-median(myktab$EndL[mySbreak:myEbreak])
    mymed_case_L[j,"startmed"]<-mystart_med
    mymed_case_L[j,"endmed"]<-myend_med
    mymed_case_L[j,"count"]<-length(mySbreak:myEbreak)
  }
}

if(length(breakpt)==0){    #if no break point found
  mymed_case_L<-matrix(0,1,3)
  colnames(mymed_case_L)<-c("startmed","endmed","count")
  mySbreak<-1
  myEbreak<-nrow(myktab)
  mystart_med<-median(myktab$StartL[mySbreak:myEbreak])
  myend_med<-median(myktab$EndL[mySbreak:myEbreak])
  mymed_case_L[1,"startmed"]<-mystart_med
  mymed_case_L[1,"endmed"]<-myend_med
  mymed_case_L[1,"count"]<-length(mySbreak:myEbreak)
}

#processing the right flank, make a table storing the coordinate for each arrow and the size
myktab<-myktab_case[order(myktab_case$StartR,decreasing=F),]
breakpt<-which(diff(myktab$StartR)>10000)

if(length(breakpt)>0){ #if at least one break point found
  mymed_case_R<-matrix(0,length(breakpt)+1,3)
  colnames(mymed_case_R)<-c("startmed","endmed","count")
  mybreak<-c(1,breakpt,nrow(myktab))
  for(j in 1:(length(mybreak)-1)){
    print(j)
    if(j==1){
      mySbreak<-1
      myEbreak<-mybreak[j+1]
    }else { #if the break value is not the first nor the last
      mySbreak<-mybreak[j]+1
      myEbreak<-mybreak[j+1]
    }
    mystart_med<-median(myktab$StartR[mySbreak:myEbreak])
    myend_med<-median(myktab$EndR[mySbreak:myEbreak])
    mymed_case_R[j,"startmed"]<-mystart_med
    mymed_case_R[j,"endmed"]<-myend_med
    mymed_case_R[j,"count"]<-length(mySbreak:myEbreak)
  }
}

if(length(breakpt)==0){    #if no break point found
  mymed_case_R<-matrix(0,1,3)
  colnames(mymed_case_R)<-c("startmed","endmed","count")
  mySbreak<-1
  myEbreak<-nrow(myktab)
  mystart_med<-median(myktab$StartR[mySbreak:myEbreak])
  myend_med<-median(myktab$EndR[mySbreak:myEbreak])
  mymed_case_R[1,"startmed"]<-mystart_med
  mymed_case_R[1,"endmed"]<-myend_med
  mymed_case_R[1,"count"]<-length(mySbreak:myEbreak)
}
mymed_ctrl_L<<-as.data.frame(mymed_ctrl_L)
mymed_ctrl_R<<-as.data.frame(mymed_ctrl_R)
mymed_case_L<<-as.data.frame(mymed_case_L)
mymed_case_R<<-as.data.frame(mymed_case_R)
#myreturn<-list(mymed_ctrl_L,mymed_ctrl_R,mymed_case_L,mymed_case_R)
#return(myreturn)

mymed_ctrl_L
mymed_ctrl_R
mymed_case_L
mymed_case_R

#convert to proportion of case and control genomes
mymed_ctrl_L$prop<-round(mymed_ctrl_L$count/nrow(myktab_ctrl)*100)
mymed_ctrl_R$prop<-round(mymed_ctrl_R$count/nrow(myktab_ctrl)*100)
mymed_case_L$prop<-round(mymed_case_L$count/nrow(myktab_case)*100)
mymed_case_R$prop<-round(mymed_case_R$count/nrow(myktab_case)*100)

###Making the plot###

#making the background plot
plot(1, type="n", xlim=c(1,5000), ylim=c(-250,250), xlab="genome position (thousands)",ylab="genome_index")

abline(h=0) #to separate the case from control kmers
text(5000,50,"case",cex=1)
text(5000,-50,"control",cex=1)


case_count<-30 #set the starting y axis level
ctrl_count<--30 #set the starting y axis level

title(main = mykmer)

#R flank (downstream) is red, L flank (upstream) is blue

#plot ctrl
for(k in 1:nrow(mymed_ctrl_L)){
  if(mymed_ctrl_L$endmed[k]<mymed_ctrl_L$startmed[k]){   #check if this flank has flipped (in the context of fwd_k)
    plot_ctrl_L_start<-mymed_ctrl_L$startmed[k]/1000+50
    plot_ctrl_L_end<-mymed_ctrl_L$endmed[k]/1000-50
  }else{
    plot_ctrl_L_start<-mymed_ctrl_L$startmed[k]/1000-50
    plot_ctrl_L_end<-mymed_ctrl_L$endmed[k]/1000+50
  }
  sizefactor<-mymed_ctrl_L$prop[k]
  x <- c(plot_ctrl_L_start, plot_ctrl_L_start, plot_ctrl_L_end)  
  #y format (top point, bottom point, central point)
  y <- c(ctrl_count-sizefactor, ctrl_count,ctrl_count-(sizefactor/2)) 
  polygon(x, y,border = NA,col=rgb(0, 0, 1,0.5))
}

ctrl_count<-ctrl_count-max(mymed_ctrl_L$prop)-30 #set the starting y axis level


for(k in 1:nrow(mymed_ctrl_R)){
  if(mymed_ctrl_R$endmed[k]<mymed_ctrl_R$startmed[k]){   #check if this flank has flipped (in the context of fwd_k)
    plot_ctrl_R_start<-mymed_ctrl_R$startmed[k]/1000+50
    plot_ctrl_R_end<-mymed_ctrl_R$endmed[k]/1000-50
  }else{
    plot_ctrl_R_start<-mymed_ctrl_R$startmed[k]/1000-50
    plot_ctrl_R_end<-mymed_ctrl_R$endmed[k]/1000+50
  }
  sizefactor<-mymed_ctrl_R$prop[k]
  x <- c(plot_ctrl_R_start, plot_ctrl_R_start, plot_ctrl_R_end)
  #y format (top point, bottom point, central point)
  y <- c(ctrl_count-sizefactor, ctrl_count,ctrl_count-(sizefactor/2))  
  polygon(x, y,border = NA,col=rgb(1, 0, 0,0.5))
}



#plot case
for(k in 1:nrow(mymed_case_L)){
  if(mymed_case_L$endmed[k]<mymed_case_L$startmed[k]){  #check if this flank has flipped (in the context of fwd_k)
    plot_case_L_start<-mymed_case_L$startmed[k]/1000+50
    plot_case_L_end<-mymed_case_L$endmed[k]/1000-50
  }else{
    plot_case_L_start<-mymed_case_L$startmed[k]/1000-50
    plot_case_L_end<-mymed_case_L$endmed[k]/1000+50
  }
  sizefactor<-mymed_case_L$prop[k]
  x <- c(plot_case_L_start, plot_case_L_start, plot_case_L_end)
  #y format (bottom point, top point, central point)
  y <- c(case_count, case_count+sizefactor,case_count+(sizefactor/2)) 
  polygon(x, y,border = NA,col=rgb(0, 0, 1,0.5))
}

case_count<-case_count+max(mymed_case_L$prop)+30 #set the starting y axis level


#case right flank arrow
for(k in 1:nrow(mymed_case_R)){
  if(mymed_case_R$endmed[k]<mymed_case_R$startmed[k]){   #check if this flank has flipped (in the context of fwd_k)
    plot_case_R_start<-mymed_case_R$startmed[k]/1000+50
    plot_case_R_end<-mymed_case_R$endmed[k]/1000-50
  }else{
    plot_case_R_start<-mymed_case_R$startmed[k]/1000-50
    plot_case_R_end<-mymed_case_R$endmed[k]/1000+50
  }
  sizefactor<-mymed_case_R$prop[k]
  x <- c(plot_case_R_start, plot_case_R_start, plot_case_R_end)  
  y <- c(case_count, case_count+sizefactor,case_count+(sizefactor/2))
  polygon(x, y,border = NA,col=rgb(1, 0, 0,0.5))
}





