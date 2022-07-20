#Usage: Rscript plot_flk_kmer_prop.R --kmer kmer93 --phen /home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv \
#--coor /home/ubuntu/Dorothy/genome_rearrangement/output/myflk_behave_pheno.txt \
#--genome.size 4000 --outdir /home/ubuntu/Dorothy/genome_rearrangement/output


#this script plot one kmer at a time
library("optparse")

option_list = list(
  make_option("--kmer", type="character", default=NULL, 
              help="kmer to plot", metavar="character"),
  make_option("--phen", type="character", default=NULL, 
              help="phenotype file, no header", metavar="character"),
  make_option("--coor", type="character", default=NULL, 
              help="kmer genome combination blast coordinates", metavar="character"),
  make_option("--genome.size", type="character", default=NULL, 
              help="size of genome in Mb", metavar="character"),
  make_option("--flk.dist", type="character", default=10000, 
              help="Maximum distance between flanks to define intact kmer", metavar="character"),
  make_option("--outdir", type="character", default=NULL,
              help="kmer to plot", metavar="character")

) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 


#setwd("/Users/dorothytam/Documents/Bath_postdoc/genomic_analyses/170_USAgenomes_analysis/111_yearGWAS//")
mypheno_file<-read.table(opt$phen,header=F)
#mypheno_file<-read.table("phenotype.txt",header=F)
#table(mypheno_file[,2])
#0  1 
#44 67

#mytable<-read.table("111yearGWAS_mystartendLR_k200.txt",header=T)
mytable<-read.table(opt$coor,header=T)
#mymerge<-merge(mytable,mypheno_file,by.x="genome",by.y="V1")
#colnames(mymerge)[7]<-"case_control"
mymerge<-mytable[,c(1,3,4,5,6,7,2)]

#hist(mytable$EndR)

#mykmer="kmer93"  #pick the kmer to plot
mykmer=opt$kmer
mykmer_tab<-mymerge[which(mymerge$kmer==mykmer),]

#making the background plot
pdf(file = paste0(opt$outdir,"/",mykmer,"_plot.pdf"),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches

x_length=(as.numeric(as.character(opt$genome.size))+1000)
#x_length=5000
plot(1, type="n", xlim=c(1,x_length), ylim=c(-350,350), xlab="genome position (thousands)")

abline(h=0) #to separate the case from control kmers
text(x_length-100,50,"case",cex=1)
text(x_length-100,-50,"control",cex=1)

title(main = paste(mykmer,"(height of arrow is the proportion of case/control genomes)",sep="\n"))

#R flank (downstream) is red, L flank (upstream) is blue

#looping through pheno and flank to plot
pheno_list<-c("ctrl","case")
flk_list<-c("L","R","in")

for(i in 1:2){
  mypheno<-pheno_list[i]
  print(mypheno)
  if(mypheno=="ctrl"){
    x=0  #for extracting pheno row from table
    startlevel=-30 #for the plotting
  }
  if(mypheno=="case"){
    x=1  #for extracting pheno row from table
    startlevel=30  #for the plotting
  }
  for(j in 1:3){
    myflk<-flk_list[j]
    print(myflk)
    if(myflk=="L"){
      mystart<-"StartL"
      myend<-"EndL"
      mycol=rgb(0, 0, 1,0.5) #blue
    }
    if(myflk=="R"){
      mystart<-"StartR"
      myend<-"EndR"
      mycol=rgb(1, 0, 0,0.5) #red
    }
    if(myflk=="in"){
      mystart<-"StartL"
      myend<-"EndR"
      mycol=rgb(0,1,0,0.5) #green
    }
    if(myflk%in%c("L","R")){  #extract the split kmer rows and specify the pheno to plot
      myktab<-mykmer_tab[which(abs(mykmer_tab$EndL-mykmer_tab$StartR)>as.numeric(as.character(opt$flk.dist)) & mykmer_tab$case_control==x),]
    }
    if(myflk=="in"){ #extract the intact kmer rows and specify the pheno to plot
      myktab<-mykmer_tab[which(abs(mykmer_tab$EndL-mykmer_tab$StartR)<as.numeric(as.character(opt$flk.dist)) & mykmer_tab$case_control==x),]
    }
    
    #make a table storing the coordinate for each arrow and the size
    print(myktab)    
    myktab<-myktab[order(myktab[[mystart]],decreasing=F),]
    breakpt<-which(diff(myktab[[mystart]])>10000)
    print(breakpt)
    
    if(length(breakpt)>0){ #if at least one break point found
      mymed<-matrix(0,length(breakpt)+1,3)
      colnames(mymed)<-c("startmed","endmed","count")
      mybreak<-c(1,breakpt,nrow(myktab))
      for(j in 1:(length(mybreak)-1)){
        #print(j)
        if(j==1){
          mySbreak<-1
          myEbreak<-mybreak[j+1]
        }else { #if the break value is not the first nor the last
          mySbreak<-mybreak[j]+1
          myEbreak<-mybreak[j+1]
        }
        mystart_med<-median(myktab[[mystart]][mySbreak:myEbreak])
        myend_med<-median(myktab[[myend]][mySbreak:myEbreak])
        mymed[j,"startmed"]<-mystart_med
        mymed[j,"endmed"]<-myend_med
        mymed[j,"count"]<-length(mySbreak:myEbreak)
      }
    }
    
    if(length(breakpt)==0){    #if no break point found
      print("break=0")     
      mymed<-matrix(0,1,3)
      colnames(mymed)<-c("startmed","endmed","count")
      mySbreak<-1
      myEbreak<-nrow(myktab)
      mystart_med<-median(myktab[[mystart]][mySbreak:myEbreak])
      myend_med<-median(myktab[[myend]][mySbreak:myEbreak])
      mymed[1,"startmed"]<-mystart_med
      mymed[1,"endmed"]<-myend_med
      mymed[1,"count"]<-length(mySbreak:myEbreak)
    }
    mymed<-as.data.frame(mymed)
    #convert to proportion of case or control genomes
    mymed$prop<-round((mymed$count/(length(which(mypheno_file[,2]==x)))*100))
    print(mymed)
    
    #plot
    for(k in 1:nrow(mymed)){
      if(mymed$endmed[k]<mymed$startmed[k]){   #check if this flank has flipped (in the context of fwd_k)
        plot_start<-mymed$startmed[k]/1000+50
        plot_end<-mymed$endmed[k]/1000-50
      }else{
        plot_start<-mymed$startmed[k]/1000-50
        plot_end<-mymed$endmed[k]/1000+50
      }
      sizefactor<-mymed$prop[k]
      
      if(mypheno=="ctrl"){
        plotx <- c(plot_start, plot_start, plot_end)  
        #y format (top point, bottom point, central point)
        ploty <- c((startlevel-sizefactor), startlevel,(startlevel-(sizefactor/2))) 
        polygon(plotx, ploty,border = NA,col=mycol)
      }
      
      if(mypheno=="case"){
        plotx <- c(plot_start, plot_start, plot_end)  
        #y format (top point, bottom point, central point)
        ploty <- c((startlevel+sizefactor), startlevel,(startlevel+(sizefactor/2))) 
        polygon(plotx, ploty,border = NA,col=mycol)
      }
    }
    
    if(mypheno=="ctrl"){
      startlevel<-(startlevel)-max(mymed$prop)-50 #set again the starting y axis level
    }
    if(mypheno=="case"){
      startlevel<-(startlevel)+max(mymed$prop)+50 #set again the starting y axis level
    }
if(myflk=="L"){
myflk_name="upstreamflk"
}
if(myflk=="R"){
myflk_name="downstreamflk"
}
if(myflk=="in"){
myflk_name="intactk"
}
    write.table(mymed,file=paste(opt$outdir,"/",mypheno,"_",myflk_name,".txt",sep=""),quote=F,col.names=T,row.names=F,sep="\t")
  }
}

#add legend
polygon(c((x_length-750),(x_length-750),(x_length-700)), c(350,330,340),border = NA,col=rgb(1, 0, 0,0.5))
text(x=(x_length-300),y=340,"kmer dowstream flank",cex = 0.6)
polygon(c((x_length-750),(x_length-750),(x_length-700)), c(320,300,310),border = NA,col=rgb(0, 0, 1,0.5))
text(x=(x_length-300),y=310,"kmer upstream flank",cex = 0.6)
polygon(c((x_length-750),(x_length-750),(x_length-700)), c(290,270,280),border = NA,col=rgb(0, 1, 0,0.5))
text(x=(x_length-450),y=280,"intact kmer",cex = 0.6)

#export the plot
dev.off()



