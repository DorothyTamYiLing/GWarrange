#input of the plotting script must be a file containing all the behaviour_summary.txt of all kmers

####################################################
setwd("/Users/dorothytam/Documents/Bath_postdoc/genomic_analyses/170_USAgenomes_analysis/")
myflk_behave_pheno<-read.table("myflk_behave_pheno.txt",header=T)
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"

#making summary information for each kmer
'%!in%' <- function(x,y)!('%in%'(x,y))

#make the final output
myall_out<-matrix(0,1,8)
colnames(myall_out)<-c("kmer","event_sum","flk_behaviour","notes","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop")

#extract the unique kmer
myk4plot<-unique(myflk_behave_pheno$kmer)


for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  #for (j in 1:10){ #open bracket for looping through each kmer
  
  mykmer<-myk4plot[j] 
  
  #select the rows referring to the kmer
  mytable<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]
  
  #making the output matrix
  myout<-matrix(0,1,8)
  colnames(myout)<-c("kmer","event_sum","flk_behaviour","notes","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop")
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
    
    #define the case and control associated "notes"
    if(mygp_ctrl_prop>0.6){   
      myout$ctrl_assos<-mygp
      myout$ctrl_assos_prop<-mygp_ctrl_prop
    }
    if(mygp_case_prop>0.6){    
      myout$case_assos<-mygp
      myout$case_assos_prop<-mygp_case_prop
    }
  }
  
  
  myout$notes<-mysum_str 
  myall_out<-rbind(myall_out,myout)
}

write.table(myall_out,file="myall_out.txt",quote=F,row.names = F,col.names = T,sep="\t")

################# below run in MAC #####################

#export  myflk_behave_pheno.txt and myall_out.txt to MAC for plotting

myall_out<-read.table("myall_out.txt",header=T)

#renaming the table
myall_out_all<-myall_out
#getting the list of kmers with event_sum!=intact_k and subsetting the corresponding rows in myflk_behave_pheno for plot
myk4plot<-myall_out[which(myall_out$event_sum!="intact_k" & myall_out$event_sum!=0),"kmer"]
myall_out<-myall_out[which(myall_out$kmer%in%myk4plot),]


myflk_behave_pheno<-read.table("myflk_behave_pheno.txt",header=T)
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"


#finding out how many kmers need to be plotting on the "case" and "control" panel
case_num<-length(which(myall_out$case_control=="1"))
ctrl_num<-length(which(myall_out$case_control=="0"))


#List of kmers for plotting
myk4plot<-unique(myall_out$kmer)

#plotting the kmers
plot(1, type="n", xlim=c(1,5000), ylim=c((-ctrl_num-(length(myk4plot)*400)),(case_num+(length(myk4plot)*400))), xlab="genome position (thousands)",ylab="genome_index")
plot(1, type="n", xlim=c(1,5000), ylim=c(-10000,10000), xlab="genome position (thousands)",ylab="genome_index")

abline(h=0) #to separate the case from control kmers
text(5000,(length(myk4plot)*10),"case",cex=0.7)
text(5000,-(length(myk4plot)*10),"control",cex=0.7)



#read in the all_behaviour_summary.txt and loop through each kmer

#kmer<-"kmer100"


#set the colour spectrum according to the number of kmers
plotcolors <- colorRampPalette(c("green","blue","red"))(length(myk4plot))
plotcolors <- colorRampPalette(c("green","blue","red"))(50)

case_count<-50
ctrl_count<--50

for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  #for (j in 1:20){  #open bracket for looping through each kmer
    
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
  points(pch=19,cex=0.5,EndL,case_count)
  }
  if(StartR>EndR){ #indicate the reverse kmer
    points(pch=19,cex=0.5,EndR,case_count)
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
    points(pch=19,cex=0.5,EndL,ctrl_count)
  }
  if(StartR>EndR){ #indicate the reverse kmer
    points(pch=19,cex=0.5,EndR,ctrl_count)
  }
  ctrl_count<-ctrl_count-1
}

    if(i==(round(nrow(mytable)/2))){
      text(10,ctrl_count,paste(myk4plot[j]," (",myctrl_assos_prop,")",sep=""),cex=0.5)
      text(10,case_count,paste(myk4plot[j]," (",mycase_assos_prop,")",sep=""),cex=0.5)
    }
    
  } #close bracket for looping through each row in the table
  
  #make the space between each kmer
  case_count=case_count+300
  ctrl_count=ctrl_count-300
  
} #close bracket for looping through each kmer

###############################################


