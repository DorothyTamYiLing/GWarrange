#(R)
#make sure the input file kmer_blastcoor.txt is in the working directory
mytable<-read.delim("kmer_blastcoor.txt", header=T,sep=" ")

#get the kmer name
kmer=mytable[1,1]
print(kmer)

#mytable<-read.table(opt$input, header=T)
#mysamplelist<-read.table(opt$list, header=T)
myphenofile<-read.table("/home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/case_control.txt", header=F)

#set the output matrix
myout<-matrix(0,0,10)
colnames(myout)<-c("case_control","genome","deletion","flk_behaviour","notes","StartL","EndL","StartR","EndR","flk_dist")

#set the output for the rows with issues
myspecial_row<-matrix(0,0,ncol(mytable))
colnames(myspecial_row)<-colnames(mytable)

#set the output for undefined behaviour of the flanking sequences
my_new_case<-matrix(0,0,ncol(mytable))
colnames(my_new_case)<-colnames(mytable)

for (i in 1:nrow(myphenofile)){
   #clear all variables
   mysample=mypheno=mysub=mynohit=StartL=EndL=StartR=EndR=mybehave=mynote="no_val"

   mysample<-as.character(myphenofile[i,1])  #get the genome name
   mypheno<-myphenofile[i,2]  #get the case/control status
   print(paste(mysample,mypheno,sep=" "))
   mysub<-mytable[which(mytable$genome==mysample),]  #extracting the rows referring to the sample, every samples should have two rows
   
   #output the sample's rows with issues (multi-hit, alignment coverage != flank length, number of SNP or gap >1) for investigation
   if (any(c("flank_needs_investigation","multi_hit","no_hit")%in%mysub$flk_presence)){   
     print("special row")
     myspecial_row<-rbind(myspecial_row,mysub)   
   }
     

   #if case of the presence of both flank, detect the flank behaviours
   if (length(mysub$flk_presence)==2 & all(mysub$flk_presence=="flank_present")){    #only do the following if both flanks are present
   
     mynohit<-"no_deletion" 
     #extract the startL, startR, endL, endR numbers
     StartL<-mysub[1,10]
     EndL<-mysub[1,11]
     StartR<-mysub[2,10]
     EndR<-mysub[2,11]
     
     #test for behaviour 1, intact kmer, case 1 forward
     if ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) < 2200)){
     mybehave<-"intact_k"
     mynote<-"fwd_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(StartR-EndL)
     } else if 
     
     #test for behaviour 1, intact kmer, case 2 reverse
     ((StartL > EndL) & (EndL > StartR) & (StartR > EndR) & ((EndL-StartR) < 2200)){
     mybehave<-"intact_k"
     mynote<-"rev_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(EndL-StartR)
     }else if
     
     #test for behaviour 2, Flank sequence move away from each other, case 1 forward 
     ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) > 2200)){
     mybehave<-"mv_away"
     mynote<-"fwd_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(StartR-EndL)
     } else if
     
     #test for behaviour 2, Flank sequence move away from each other, case 2 reverse
     ((StartL > EndL) & (EndL > StartR) & (StartR > EndR) & ((EndL-StartR) > 2200)){
     mybehave<-"mv_away"
     mynote<-"rev_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(EndL-StartR)
     } else if
     
     #test for behaviour 3, Left flank and right flank swap position
     ((StartR < EndR) & (StartL < EndL)  & (EndR < StartL)){
     mybehave<-"swp_flk"
     mynote<-"NA"
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(StartL-EndR)
     } else if
     
     #test for behaviour 4, flank sequence move away from each other and flipped, case 1, right flank move and flipped
     ((StartL < EndL) & (EndL < EndR) & (EndR < StartR)){
     mybehave<-"mv&flp"
     mynote<-"mv&flp"   #not enough information to tell whether it is the right or left flank that is moved and flipped, need to see is the intact_k of the case/control counterpart is forward or reversed
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(EndR-EndL)
     } else if
     
     #test for behaviour 4, flank sequence move away from each other and flipped, case 1, left flank move and flipped
     ((EndL < StartL) & (StartL < StartR) & (StartR < EndR)){
     mybehave<-"mv&flp"
     mynote<-"mv&flp"    #not enough information to tell whether it is the right or left flank that is moved and flipped, need to see is the intact_k of the case/control counterpart is forward or reversed
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(StartR-StartL)
     } else if
     
     #test for behaviour 5, flank sequence swap position and flipped, case 1, left flank move and flipped
     ((StartR < EndR) & (EndR < EndL) & (EndL < StartL)){
     mybehave<-"swp&flp"
     mynote<-"swp&flp"   #not enough information to tell whether it is the right or left flank that is swapped and flipped, need to see is the intact_k of the case/control counterpart is forward or reversed
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(EndL-EndR)
     } else if
     
     #test for behaviour 5, flank sequence swap position and flipped, case 2, right flank move and flipped
     ((EndR < StartR) & (StartR < StartL) & (StartL < EndL)){
     mybehave<-"swp&flp"
     mynote<-"swp&flp"    #not enough information to tell whether it is the right or left flank that is swapped and flipped, need to see is the intact_k of the case/control counterpart is forward or reversed
     StartL_out<-StartL
     EndL_out<-EndL
     StartR_out<-StartR
     EndR_out<-EndR
     flk_dist<-(StartL-StartR)
     } else {
     
     #output undefined behaviour of the flanking sequences
     my_new_case<-rbind(my_new_case, mysub)
     }
     myoutrow<-c(mypheno,mysample,mynohit,mybehave,mynote,StartL_out,EndL_out,StartR_out,EndR_out,flk_dist)
   myout<-rbind(myout,myoutrow)  #append the info for output
   }   #closing bracket of the "all flank_present" 
}   #closing bracket for the samplelist for-loop

write.table(myout,file="kmer_behaviour_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")
if (nrow(myspecial_row)>0){
write.table(myspecial_row,file=paste(kmer,"flank_with_issues.txt",sep="_"),quote=F,row.names = F,col.names = T,sep="\t")
}
if (nrow(my_new_case)>0){
write.table(my_new_case,file=paste(kmer,"new_flk_behaviour.txt",sep="_"),quote=F,row.names = F,col.names = T,sep="\t")
}

####################################################

#only do the following if the number of rows in myout>1
if(nrow(myout)>0){ 


'%!in%' <- function(x,y)!('%in%'(x,y))

mytable<-read.table("kmer_behaviour_summary.txt", header=T)

#making the output matrix
myout<-matrix(0,1,17)
colnames(myout)<-c("deletion","flk_behaviour","notes","StartL","EndL","StartR","EndR","flk_dist", "summary_deletion","flk_summary","casectrl","leftright","StartL_med","EndL_med","StartR_med","EndR_med","prop_comp")
myout<-as.data.frame(myout)


#get the total number of cases and controls
ctrl_count<-length(which(mytable$case_control=="0"))
case_count<-length(which(mytable$case_control=="1"))

#count the proportion of cases and controls in deletion column
mytable$deletion<-as.character(mytable$deletion) 
mycat<-unique(mytable$deletion)
mysum_str<-""
for (i in 1:length(mycat)){
mygp<-as.character(mycat[i])  #extract the group

#count number and proportion of gp in ctrl genomes
mygp_ctrl<-length(which(mytable$deletion==mygp & mytable$case_control=="0"))
mygp_ctrl_prop<-round(mygp_ctrl/ctrl_count,2)
myctrl_str<-paste(mygp_ctrl,"(",mygp_ctrl_prop,")",sep="")

#count number and proportion of gp in case genomes
mygp_case<-length(which(mytable$deletion==mygp & mytable$case_control=="1"))
mygp_case_prop<-round(mygp_case/case_count,2)
mycase_str<-paste(mygp_case,"(",mygp_case_prop,")",sep="")

mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make a string that show deletion status:control_count:case_count
mysum_str<-paste(mysum_str,mysum,sep=" ")
}

myout$deletion<-mysum_str

#making summary for deletion 
if(all(c(mygp_ctrl_prop,mygp_case_prop)==1)){
myout$summary_deletion<-"no_deletion"
}else{
myout$summary_deletion<-"with_deletion"
}

#get the summary statistics of the StartL in case and control genomes
ctrl_StartL_stat<-round(summary(mytable[which(mytable$case_control=="0"),6]),0)
ctrl_StartL_stat_str_1<-paste(ctrl_StartL_stat, collapse = ' ')  #make summary stat string
ctrl_StartL_SD<-round(sd(mytable[which(mytable$case_control=="0"),6]),0)
ctrl_StartL_stat_str<-paste(ctrl_StartL_stat_str_1, ctrl_StartL_SD,sep=" ") #make new string including SD
ctrl_StartL_med<-round(median(mytable[which(mytable$case_control=="0"),6]),0)

case_StartL_stat<-round(summary(mytable[which(mytable$case_control=="1"),6]),0)
case_StartL_stat_str_1<-paste(case_StartL_stat, collapse = ' ') #make summary stat string
case_StartL_SD<-round(sd(mytable[which(mytable$case_control=="1"),6]),0)
case_StartL_stat_str<-paste(case_StartL_stat_str_1, case_StartL_SD,sep=" ") #make new string including SD
case_StartL_med<-round(median(mytable[which(mytable$case_control=="1"),6]),0)

myout$StartL<-paste(case_StartL_stat_str,ctrl_StartL_stat_str,sep=" | ")
myout$StartL_med<-paste(case_StartL_med,ctrl_StartL_med,sep=" | ")

#get the summary statistics of the EndL in case and control genomes
ctrl_EndL_stat<-round(summary(mytable[which(mytable$case_control=="0"),7]),0)
ctrl_EndL_stat_str_1<-paste(ctrl_EndL_stat, collapse = ' ')  #make summary stat string
ctrl_EndL_SD<-round(sd(mytable[which(mytable$case_control=="0"),7]),0)
ctrl_EndL_stat_str<-paste(ctrl_EndL_stat_str_1, ctrl_EndL_SD,sep=" ") #make new string including SD
ctrl_EndL_med<-round(median(mytable[which(mytable$case_control=="0"),7]),0)

case_EndL_stat<-round(summary(mytable[which(mytable$case_control=="1"),7]),0)
case_EndL_stat_str_1<-paste(case_EndL_stat, collapse = ' ') #make summary stat string
case_EndL_SD<-round(sd(mytable[which(mytable$case_control=="1"),7]),0)
case_EndL_stat_str<-paste(case_EndL_stat_str_1, case_EndL_SD,sep=" ") #make new string including SD
case_EndL_med<-round(median(mytable[which(mytable$case_control=="1"),7]),0)

myout$EndL<-paste(case_EndL_stat_str,ctrl_EndL_stat_str,sep=" | ")
myout$EndL_med<-paste(case_EndL_med,ctrl_EndL_med,sep=" | ")

#get the summary statistics of the StartR in case and control genomes
ctrl_StartR_stat<-round(summary(mytable[which(mytable$case_control=="0"),8]),0)
ctrl_StartR_stat_str_1<-paste(ctrl_StartR_stat, collapse = ' ')  #make summary stat string
ctrl_StartR_SD<-round(sd(mytable[which(mytable$case_control=="0"),8]),0)
ctrl_StartR_stat_str<-paste(ctrl_StartR_stat_str_1, ctrl_StartR_SD,sep=" ") #make new string including SD
ctrl_StartR_med<-round(median(mytable[which(mytable$case_control=="0"),8]),0)

case_StartR_stat<-round(summary(mytable[which(mytable$case_control=="1"),8]),0)
case_StartR_stat_str_1<-paste(case_StartR_stat, collapse = ' ') #make summary stat string
case_StartR_SD<-round(sd(mytable[which(mytable$case_control=="1"),8]),0)
case_StartR_stat_str<-paste(case_StartR_stat_str_1, case_StartR_SD,sep=" ") #make new string including SD
case_StartR_med<-round(median(mytable[which(mytable$case_control=="1"),8]),0)

myout$StartR<-paste(case_StartR_stat_str,ctrl_StartR_stat_str,sep=" | ")
myout$StartR_med<-paste(case_StartR_med,ctrl_StartR_med,sep=" | ")

#get the summary statistics of the EndR in case and control genomes
ctrl_EndR_stat<-round(summary(mytable[which(mytable$case_control=="0"),9]),0)
ctrl_EndR_stat_str_1<-paste(ctrl_EndR_stat, collapse = ' ')  #make summary stat string
ctrl_EndR_SD<-round(sd(mytable[which(mytable$case_control=="0"),9]),0)
ctrl_EndR_stat_str<-paste(ctrl_EndR_stat_str_1, ctrl_EndR_SD,sep=" ") #make new string including SD
ctrl_EndR_med<-round(median(mytable[which(mytable$case_control=="0"),9]),0)

case_EndR_stat<-round(summary(mytable[which(mytable$case_control=="1"),9]),0)
case_EndR_stat_str_1<-paste(case_EndR_stat, collapse = ' ') #make summary stat string
case_EndR_SD<-round(sd(mytable[which(mytable$case_control=="1"),9]),0)
case_EndR_stat_str<-paste(case_EndR_stat_str_1, case_EndR_SD,sep=" ") #make new string including SD
case_EndR_med<-round(median(mytable[which(mytable$case_control=="1"),9]),0)

myout$EndR<-paste(case_EndR_stat_str,ctrl_EndR_stat_str,sep=" | ")
myout$EndR_med<-paste(case_EndR_med,ctrl_EndR_med,sep=" | ")

#get the summary statistics of the flk_dist in case and control genomes
ctrl_flk_dist_stat<-round(summary(mytable[which(mytable$case_control=="0"),10]),0)
ctrl_flk_dist_stat_str<-paste(ctrl_flk_dist_stat, collapse = ' ')
case_flk_dist_stat<-round(summary(mytable[which(mytable$case_control=="1"),10]),0)
case_flk_dist_stat_str<-paste(case_flk_dist_stat, collapse = ' ')
myout$flk_dist<-paste(case_flk_dist_stat_str,ctrl_flk_dist_stat_str,sep=" | ")

#count the proportion of cases and controls in flk_behaviour column
mytable$flk_behaviour<-as.character(mytable$flk_behaviour) 
mycat<-unique(mytable$flk_behaviour)
mysum_str<-""

#store the groups in a list
mygplist<-c()

for (i in 1:length(mycat)){
mygp<-as.character(mycat[i])  #extract the group
mygplist<-c(mygplist,mygp)

#count number and proportion of gp in ctrl genomes
mygp_ctrl<-length(which(mytable$flk_behaviour==mygp & mytable$case_control=="0"))
mygp_ctrl_prop_FB<-round(mygp_ctrl/ctrl_count,2)
myctrl_str<-paste(mygp_ctrl,"(",mygp_ctrl_prop_FB,")",sep="")

#count number and proportion of gp in case genomes
mygp_case<-length(which(mytable$flk_behaviour==mygp & mytable$case_control=="1"))
mygp_case_prop_FB<-round(mygp_case/case_count,2)
mycase_str<-paste(mygp_case,"(",mygp_case_prop_FB,")",sep="")

mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make a string that show deletion status:control_count:case_count
mysum_str<-paste(mysum_str,mysum,sep=" ")
}

myout$flk_behaviour<-mysum_str


#count the proportion of cases and controls in notes column
mytable$notes<-as.character(mytable$notes) 
mycat<-unique(mytable$notes)
mysum_str<-""
for (i in 1:length(mycat)){
mygp<-as.character(mycat[i])  #extract the group

#count number and proportion of gp in ctrl genomes
mygp_ctrl<-length(which(mytable$notes==mygp & mytable$case_control=="0"))
mygp_ctrl_prop_note<-round(mygp_ctrl/ctrl_count,2)
myctrl_str<-paste(mygp_ctrl,"(",mygp_ctrl_prop_note,")",sep="")

#count number and proportion of gp in case genomes
mygp_case<-length(which(mytable$notes==mygp & mytable$case_control=="1"))
mygp_case_prop_note<-round(mygp_case/case_count,2)
mycase_str<-paste(mygp_case,"(",mygp_case_prop_note,")",sep="")

mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make a string that show deletion status:control_count:case_count
mysum_str<-paste(mysum_str,mysum,sep=" ")
}

myout$notes<-mysum_str   

#from now on mysum_str refers to notes column only

#getting flank behaviour and note summary for the kmer (summary of flank behaviour, case/control, left/right reversed)
#for intact_k only, the 
if(all(c("mv_away","swp_flk","mv&flp","swp&flp")%!in%mygplist)){
myout$flk_summary<-"intact_k only"
if(grepl("rev_k",mysum_str,fixed = TRUE)){ # if there is reversed kmer
mysplit<-unlist(strsplit(mysum_str," "))
myrevk<-mysplit[which(grepl("rev_k",mysplit,fixed = TRUE))]  #extract the rev kmer section
myrevk<-gsub("(", ":", myrevk,fixed = TRUE) #replace "("
myrevk<-gsub(")", "", myrevk,fixed = TRUE)   #replace ")"
myrevk<-unlist(strsplit(myrevk,":"))  #get the numbers 
if(myrevk[3]=="1" & myrevk[5]=="0"){
myrev_gp<-"case"     #deciding if case and control are 
}
if(myrevk[5]=="1" & myrevk[3]=="0"){
myrev_gp<-"control"
}
if(myrevk[5]=="1" & myrevk[3]=="1"){
myrev_gp<-"case and control"
}
myout$casectrl<-myrev_gp
}else{
myout$casectrl<-"no_rev_k"  #if there is not any reverse kmer
}
myout$leftright<-NA
}

#for mv_away
if("mv_away"%in%mygplist){
myout$flk_summary<-"mv_away"
mysplit<-unlist(strsplit(mysum_str," "))
myflp<-mysplit[which(grepl("mv_away",mysplit,fixed = TRUE))]  #extract the rev kmer section
myflp<-gsub("(", ":", myflp,fixed = TRUE) #replace "("
myflp<-gsub(")", "", myflp,fixed = TRUE)   #replace ")"
myflp<-unlist(strsplit(myflp,":"))  #get the numbers 
if(myflp[3]=="1"){
myflp_gp<-"case"     #deciding if case and control are 
}
if(myflp[5]=="1"){
myflp_gp<-"control"
}
myout$casectrl<-myflp_gp
Ldiff<-abs(ctrl_StartL_med-case_StartL_med)
Rdiff<-abs(ctrl_StartR_med-case_StartR_med)
if(Ldiff>Rdiff){
myout$leftright<-"Lflk"
}else{
myout$leftright<-"Rflk"
}
}

#for swp_flk
if("swp_flk"%in%mygplist){
myout$flk_summary<-"swp_flk"
mysplit<-unlist(strsplit(mysum_str," "))
myflp<-mysplit[which(grepl("swp_flk",mysplit,fixed = TRUE))]  #extract the rev kmer section
myflp<-gsub("(", ":", myflp,fixed = TRUE) #replace "("
myflp<-gsub(")", "", myflp,fixed = TRUE)   #replace ")"
myflp<-unlist(strsplit(myflp,":"))  #get the numbers 
if(myflp[3]=="1"){
myflp_gp<-"case"     #deciding if case and control are 
}
if(myflp[5]=="1"){
myflp_gp<-"control"
}
myout$casectrl<-myflp_gp
Ldiff<-abs(ctrl_StartL_med-case_StartL_med)
Rdiff<-abs(ctrl_StartR_med-case_StartR_med)
if(Ldiff>Rdiff){
myout$leftright<-"Lflk"
}else{
myout$leftright<-"Rflk"
}
}


#for mv&flp 
if("mv&flp"%in%mygplist){
myout$flk_summary<-"mv&flp"
mysplit<-unlist(strsplit(mysum_str," "))
myflp<-mysplit[which(grepl("mv&flp",mysplit,fixed = TRUE))]  #extract the rev kmer section
myflp<-gsub("(", ":", myflp,fixed = TRUE) #replace "("
myflp<-gsub(")", "", myflp,fixed = TRUE)   #replace ")"
myflp<-unlist(strsplit(myflp,":"))  #get the numbers 
if(myflp[3]=="1"){
myflp_gp<-"case"     #deciding if case and control are 
}
if(myflp[5]=="1"){
myflp_gp<-"control"
}
myout$casectrl<-myflp_gp
Ldiff<-abs(ctrl_StartL_med-case_StartL_med)
Rdiff<-abs(ctrl_StartR_med-case_StartR_med)
if(Ldiff>Rdiff){
myout$leftright<-"Lflk"
}else{
myout$leftright<-"Rflk"
}
}


#for swp&flp 
if("swp&flp"%in%mygplist){
myout$flk_summary<-"swp&flp"
mysplit<-unlist(strsplit(mysum_str," "))
myflp<-mysplit[which(grepl("swp&flp",mysplit,fixed = TRUE))]  #extract the rev kmer section
myflp<-gsub("(", ":", myflp,fixed = TRUE) #replace "("
myflp<-gsub(")", "", myflp,fixed = TRUE)   #replace ")"
myflp<-unlist(strsplit(myflp,":"))  #get the numbers 
if(myflp[3]=="1"){
myflp_gp<-"case"     #deciding if case and control are 
}
if(myflp[5]=="1"){
myflp_gp<-"control"
}
myout$casectrl<-myflp_gp
Ldiff<-abs(ctrl_StartL_med-case_StartL_med)
Rdiff<-abs(ctrl_StartR_med-case_StartR_med)
if(Ldiff>Rdiff){
myout$leftright<-"Lflk"
}else{
myout$leftright<-"Rflk"
}
}

#detecting if all the flk_behaviour proportion and notes proportion values are either "0" or "1", if any proportion detected were not 0/1, my_prop_count=my_prop_count+1
if (all(c(mygp_ctrl_prop_FB,mygp_case_prop_FB,mygp_ctrl_prop_note,mygp_case_prop_note)%in%c(0,1))){
myout$prop_comp<-"complete_prop"
}else{
myout$prop_comp<-"some_dissociation"
}

if(nrow(myout)>0){  #only output if there is information
write.table(myout,file=paste(kmer,"oneline_summary.txt",sep="_"),quote=F,col.names=T,row.names=F,sep="\t")
}
#myout_t<-t(myout)
#write.table(myout_t,file="table_summary.txt",quote=F,col.names=F,row.names=T,sep="\t")

} #closing of the nrow(myout)>0 if statement
