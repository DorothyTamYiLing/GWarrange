#input of the plotting script must be a file containing all the behaviour_summary.txt of all kmers

####################################################
setwd("/Users/dorothytam/Desktop/")
myflk_behave_pheno<-read.table("myflk_behave_pheno.txt",header=T)
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"

#remove the rows with flank behaviour == NA
myflk_behave_pheno<-myflk_behave_pheno[which(!is.na(myflk_behave_pheno$flk_behaviour)),]

#make the final output
myall_out<-matrix(0,1,14)
colnames(myall_out)<-c("kmer","event_sum","flk_behaviour","notes","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop","case_assos_gp_Lflk_sumstat","case_assos_gp_Rflk_sumstat","ctrl_assos_gp_Lflk_sumstat","ctrl_assos_gp_Rflk_sumstat","case_assos_gp_flkdis_sumstat","ctrl_assos_gp_flkdis_sumstat")

#extract the unique kmer
myk4plot<-unique(myflk_behave_pheno$kmer)


for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  #for (j in 1:10){ #open bracket for looping through each kmer
  
  mykmer<-myk4plot[j] 
  
  #select the rows referring to the kmer
  mytable<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]
  
  #making the output matrix
myout<-matrix(0,1,14)
colnames(myout)<-c("kmer","event_sum","flk_behaviour","notes","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop","case_assos_gp_Lflk_sumstat","case_assos_gp_Rflk_sumstat","ctrl_assos_gp_Lflk_sumstat","ctrl_assos_gp_Rflk_sumstat","case_assos_gp_flkdis_sumstat","ctrl_assos_gp_flkdis_sumstat")
myout<-as.data.frame(myout)
myout$kmer<-mykmer
  
  #get the total number of cases and controls
  ctrl_count<-length(which(mytable$case_control=="0"))
  case_count<-length(which(mytable$case_control=="1"))
  
  #count the proportion of cases and controls in flk_behaviour column
  mytable$flk_behaviour<-as.character(mytable$flk_behaviour) 
  mycat<-unique(mytable$flk_behaviour)
  
  myout$event_sum<-paste(mycat,collapse=":")
  
  mysum_str<-""
  
  for (i in 1:length(mycat)){
    mygp<-as.character(mycat[i])  #extract the behave group name
    
    #count number and proportion of gp in ctrl genomes
    mygp_ctrl<-length(which(mytable$flk_behaviour==mygp & mytable$case_control=="0"))
    mygp_ctrl_prop<-round(mygp_ctrl/ctrl_count,2)
    myctrl_str<-paste(mygp_ctrl,"(",mygp_ctrl_prop,")",sep="")
    
    #count number and proportion of gp in case genomes
    mygp_case<-length(which(mytable$flk_behaviour==mygp & mytable$case_control=="1"))
    mygp_case_prop<-round(mygp_case/case_count,2)
    mycase_str<-paste(mygp_case,"(",mygp_case_prop,")",sep="")
    
    mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make a string that show deletion status:control_count:case_count
    mysum_str<-paste(mysum_str,mysum,sep=" ")
  }
  
  myout$flk_behaviour<-mysum_str
  
  #clear the variables
  mycat=mysum_str=mygp=mygp_ctrl=mygp_ctrl_prop=myctrl_str=mygp_case=mygp_case_prop=mycase_str=mysum=mysum_str=NA
  
  #count the proportion of cases and controls in notes column
  mytable$notes<-as.character(mytable$notes) 
  mycat<-unique(mytable$notes)
  mysum_str<-""
  
  for (i in 1:length(mycat)){
    mygp<-as.character(mycat[i])  #extract the notes group name
    
    #count number and proportion of gp in ctrl genomes
    mygp_ctrl<-length(which(mytable$notes==mygp & mytable$case_control=="0"))
    mygp_ctrl_prop<-round(mygp_ctrl/ctrl_count,2)
    myctrl_str<-paste(mygp_ctrl,"(",mygp_ctrl_prop,")",sep="") 
    
    #count number and proportion of gp in case genomes
    mygp_case<-length(which(mytable$notes==mygp & mytable$case_control=="1"))
    mygp_case_prop<-round(mygp_case/case_count,2)
    mycase_str<-paste(mygp_case,"(",mygp_case_prop,")",sep="")
    
    mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make a string that show deletion status:control_count:case_count
    mysum_str<-paste(mysum_str,mysum,sep=" ")
    
    myout$notes<-mysum_str
    
    #define the case and control associated "notes"
    if(mygp_ctrl_prop>0.6){  
      myctrl_assos_gp<-mygp
      myout$ctrl_assos<-mygp
      myout$ctrl_assos_prop<-mygp_ctrl_prop
    }
    if(mygp_case_prop>0.6){ 
      mycase_assos_gp<-mygp   
      myout$case_assos<-mygp
      myout$case_assos_prop<-mygp_case_prop
    }
    }
    
    #after determining the case and control associated "notes"....
    
#get the coordinate summary statistics (StartL and endR) of the myctrl_assos_gp in ctrl genomes according to note groups 
myStartLstat<-paste(round(summary(mytable[which(mytable$notes==myctrl_assos_gp & mytable$case_control=="0"),"StartL"]),0),collapse=" ")
myStartLSD<-round(sd(mytable[which(mytable$notes==myctrl_assos_gp & mytable$case_control=="0"),"StartL"]),0)
myctrl_Lflk_sumstat_str<-paste(myStartLstat,myStartLSD,sep=" ")

myStartRstat<-paste(round(summary(mytable[which(mytable$notes==myctrl_assos_gp & mytable$case_control=="0"),"StartR"]),0),collapse=" ")
myStartRSD<-round(sd(mytable[which(mytable$notes==myctrl_assos_gp & mytable$case_control=="0"),"StartR"]),0)
myctrl_Rflk_sumstat_str<-paste(myStartRstat,myStartRSD,sep=" ")

myflkdiststat_str<-paste(round(summary(mytable[which(mytable$notes==myctrl_assos_gp & mytable$case_control=="0"),"flk_dist"]),0),collapse=" ")

myout$ctrl_assos_gp_Lflk_sumstat<-myctrl_Lflk_sumstat_str
myout$ctrl_assos_gp_Rflk_sumstat<-myctrl_Rflk_sumstat_str
myout$ctrl_assos_gp_flkdis_sumstat<-myflkdiststat_str

#clear the variables
myStartLstat=myStartLSD=myStartRstat=myStartRSD=myflkdiststat=NA

#get the coordinate summary statistics (StartL and endR) of the kmer in case genomes according to note groups
myStartLstat<-paste(round(summary(mytable[which(mytable$notes==mycase_assos_gp & mytable$case_control=="1"),"StartL"]),0),collapse=" ")
myStartLSD<-round(sd(mytable[which(mytable$notes==mycase_assos_gp & mytable$case_control=="1"),"StartL"]),0)
mycase_Lflk_sumstat_str<-paste(myStartLstat,myStartLSD,sep=" ")

myStartRstat<-paste(round(summary(mytable[which(mytable$notes==mycase_assos_gp & mytable$case_control=="1"),"StartR"]),0),collapse=" ")
myStartRSD<-round(sd(mytable[which(mytable$notes==mycase_assos_gp & mytable$case_control=="1"),"StartR"]),0)
mycase_Rflk_sumstat_str<-paste(myStartRstat,myStartRSD,sep=" ")

myflkdiststat_str<-paste(round(summary(mytable[which(mytable$notes==mycase_assos_gp & mytable$case_control=="1"),"flk_dist"]),0),collapse=" ")

myout$case_assos_gp_Lflk_sumstat<-mycase_Lflk_sumstat_str
myout$case_assos_gp_Rflk_sumstat<-mycase_Rflk_sumstat_str
myout$case_assos_gp_flkdis_sumstat<-myflkdiststat_str
  
  myall_out<-rbind(myall_out,myout)
}

write.table(myall_out,file="myall_out.txt",quote=F,row.names = F,col.names = T,sep="\t")


################# below run in MAC #####################

#plotting based on myflk_behave_pheno.txt only

####plotting#####

myk4plot<-unique(myflk_behave_pheno[which(myflk_behave_pheno$flk_behaviour!="intact_k"),"kmer"])

mykmer<-"kmer10"

#extract the table of the kmer
myktab<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]

abline(h=0) #to separate the case from control kmers
text(5000,(length(myk4plot)*10),"case",cex=0.7)
text(5000,-(length(myk4plot)*10),"control",cex=0.7)

case_count<-100 #set the starting y axis level
ctrl_count<--100 #set the starting y axis level

make_arrow_coor("kmer10")

mykorien<-unique(myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),"myk_orien"])
if(mykorien=="fwd_k"){
  mymed_ctrl_L$startmed<-mymed_ctrl_L$startmed/1000-50
  mymed_ctrl_R$endmed<-mymed_ctrl_R$endmed/1000+50
  mymed_case_L$startmed<-mymed_case_L$startmed/1000-50
  mymed_case_R$endmed<-mymed_case_R$endmed/1000+50
}

if(mykorien=="rev_k"){
  mymed_ctrl_L$startmed<-mymed_ctrl_L$startmed/1000+50
  mymed_ctrl_R$endmed<-mymed_ctrl_R$endmed/1000-50
  mymed_case_L$startmed<-mymed_case_L$startmed/1000+50
  mymed_case_R$endmed<-mymed_case_R$endmed/1000-50
}

#plotting the finalised arrows


arrow_coor_casectrl<-function("0","ctrl"){ #function with kmer function , making mymed_ctrl_L and mymed_ctrl_R
}
arrow_coor_casectrl<-function("1","case"){ #function with kmer function,  making mymed_case_L and mymed_case_R
}
#assign("myktab_ctr",myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer & myflk_behave_pheno$case_control=="0"),])
  

#function for making mymed_ctrl_L, mymed_ctrl_R, mymed_case_L, mymed_case_R for each kmer
make_arrow_coor<-function(mykmer){
  #extract the table of the kmer and ctrl
  myktab_ctrl<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer & myflk_behave_pheno$case_control=="0"),]
  
  
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
}




#### old code #####

#making a function
'%!in%' <- function(x,y)!('%in%'(x,y))

setwd("/Users/dorothytam/Desktop/")

myflk_behave_pheno<-read.table("myflk_behave_pheno.txt",header=T)
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"

#remove the rows with flank behaviour == NA
myflk_behave_pheno<-myflk_behave_pheno[which(!is.na(myflk_behave_pheno$flk_behaviour)),]


#select the kmers that are with at least one row of flank_behaviour!=intact_k
myk4plot<-unique(myflk_behave_pheno[which(myflk_behave_pheno$flk_behaviour!="intact_k"),"kmer"])


#case panel y axis
control_y<--(length(myk4plot)*control_num+length(myk4plot)*300)
case_y<-length(myk4plot)*case_num+length(myk4plot)*300


#plotting the kmers
#plot(1, type="n", xlim=c(1,5000), ylim=c(control_y,case_y), xlab="genome position (thousands)",ylab="genome_index")
plot(1, type="n", xlim=c(1,5000), ylim=c(-50000,50000), xlab="genome position (thousands)",ylab="genome_index")

abline(h=0) #to separate the case from control kmers
text(5000,(length(myk4plot)*10),"case",cex=0.7)
text(5000,-(length(myk4plot)*10),"control",cex=0.7)

#set the colour spectrum according to the number of kmers
#plotcolors <- colorRampPalette(c("green","blue","red"))(length(myk4plot))
plotcolors <- colorRampPalette(c("green","blue","red"))(30)

case_count<-100 #set the starting y axis level
ctrl_count<--100 #set the starting y axis level

#for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  for (j in 1:30){  #open bracket for looping through each kmer
  
  #kmer<-myk4plot[j]
  mykmer<-"kmer10"
  
  #extract the table of the kmer
  myktab<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]
  
  #processing the left flanks, then right flanks
  myktab<-myktab[order(myktab$StartL,decreasing=F),]
  
  
  #extract the case and control associated kmer "notes" information from myall_out
  mycase_assos<-myall_out[which(myall_out$kmer==myk4plot[j]),"case_assos"]
  myctrl_assos<-myall_out[which(myall_out$kmer==myk4plot[j]),"ctrl_assos"]
  
  mycase_assos_prop<-myall_out[which(myall_out$kmer==myk4plot[j]),"case_assos_prop"]
  myctrl_assos_prop<-myall_out[which(myall_out$kmer==myk4plot[j]),"ctrl_assos_prop"]
  
  
  #select the rows referring to the kmer
  mytable<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==myk4plot[j]),]
  
  #plotting each kmer mapped on each genome
  for (i in 1:nrow(mytable)){  #looping through each row in the table
    print(i)
    StartL<-as.numeric(mytable[i,"StartL"])/1000
    EndL<-as.numeric(mytable[i,"EndL"])/1000
    StartR<-as.numeric(mytable[i,"StartR"])/1000
    EndR<-as.numeric(mytable[i,"EndR"])/1000
    
    #when the flanks are apart
    if(mytable[i,"flk_behaviour"]!="intact_k"){ 
      if(StartL>EndL){
        StartL=StartL+25
        EndL<-EndL-25
      }else{
        StartL=StartL-25
        EndL<-EndL+25
      }
      if(StartR>EndR){
        StartR=StartR+25
        EndR<-EndR-25
      }else{
        StartR=StartR-25
        EndR<-EndR+25
      }
    }
    
    #when the flanks are together, i.e. intact k, and k is mapped forward
    if(mytable[i,"flk_behaviour"]=="intact_k"){ 
      if(EndR>StartL){
        StartL=StartL-50
        EndR<-EndR+50
      }
      if(StartL>EndR){
        StartL=StartL+50
        EndR<-EndR-50
      }
    }
    
    #decide if its case or control
    
    #if it is case
    if(mytable[i,"case_control"]=="1"){ 
      #if the kmer behaviour correspond to the case associated event
      if(mytable[i,"notes"]==mycase_assos){
        arrows(x0=StartL,x1=EndL, y0=case_count, y1=case_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col=plotcolors[j], add=T)  #plotting the left flank
        arrows(x0=StartR,x1=EndR, y0=case_count, y1=case_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col=plotcolors[j], add=T)  #plotting the right flank
      }
      #if the kmer behaviour correspond to the control associated event
      if(mytable[i,"notes"]==myctrl_assos){
        arrows(x0=StartL,x1=EndL, y0=case_count, y1=case_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col="grey", add=T)  #plotting the left flank
        arrows(x0=StartR,x1=EndR, y0=case_count, y1=case_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col="grey", add=T)  #plotting the right flank
      }
      if(StartL>EndL){ #indicate the reverse kmer
        #points(pch=19,cex=0.5,EndL,case_count)
      }
      if(StartR>EndR){ #indicate the reverse kmer
        #points(pch=19,cex=0.5,EndR,case_count)
      }
      case_count<-case_count+1
    }
    
    #if it is control
    if(mytable[i,"case_control"]=="0"){ 
      #if the kmer behaviour correspond to the control associated event
      if(mytable[i,"notes"]==myctrl_assos ){
        arrows(x0=StartL,x1=EndL, y0=ctrl_count, y1=ctrl_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col=plotcolors[j], add=T)  #plotting the left flank
        arrows(x0=StartR,x1=EndR, y0=ctrl_count, y1=ctrl_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col=plotcolors[j], add=T)  #plotting the right flank
      }
      #if the kmer behaviour correspond to the case associated event
      if(mytable[i,"notes"]==mycase_assos){
        arrows(x0=StartL,x1=EndL, y0=ctrl_count, y1=ctrl_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col="grey", add=T)  #plotting the left flank
        arrows(x0=StartR,x1=EndR, y0=ctrl_count, y1=ctrl_count, length =0.05, lwd=1, angle = 20,arr.width=0.1, code = 2,col="grey", add=T)  #plotting the right flank
      }
      if(StartL>EndL){ #indicate the reverse kmer
        #points(pch=19,cex=0.5,EndL,ctrl_count)
      }
      if(StartR>EndR){ #indicate the reverse kmer
        #points(pch=19,cex=0.5,EndR,ctrl_count)
      }
      ctrl_count<-ctrl_count-1
    }
    
    if(i==(round(nrow(mytable)/2))){
      text(10,ctrl_count,paste(myk4plot[j]," (",myctrl_assos_prop,")",sep=""),cex=0.3)
      text(10,case_count,paste(myk4plot[j]," (",mycase_assos_prop,")",sep=""),cex=0.3)
    }
    
  } #close bracket for looping through each row in the table
  
  #make the space between each kmer
  case_count=case_count+1500
  ctrl_count=ctrl_count-1500
  
} #close bracket for looping through each kmer

###############################################



