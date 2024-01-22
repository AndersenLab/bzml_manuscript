library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ape)
library('Biostrings')
library(grid)
library(gridExtra)
library(cowplot)
library(showtext)
showtext_auto()
showtext_opts(dpi=600,)

plotNormGeneFeatures <- function(sub_gff,cust) {
  
  GFFbyL2 <- sub_gff %>% 
    dplyr::arrange(desc(spp),L2) %>%
    dplyr::group_by(spp,L2) %>%
    dplyr::mutate(ngroup=cur_group_id()) %>%
    dplyr::ungroup()
  
  if (cust == T) {
    groups <- data.frame(mAlias=c("ben-1","avr-14","avr-15","glc-1"),ngroup=c(4,3,2,1))
    GFFbyL2 <- GFFbyL2 %>%
      dplyr::select(-ngroup) %>%
      dplyr::left_join(groups,by="mAlias")
  }
  
  neg_str <- GFFbyL2 %>% dplyr::filter(strand=="-")
  neg_L3 <- neg_str %>% 
    dplyr::filter(type=="exon" | type =="intron" | type =="CDS" | type=="DELETION") %>%
    dplyr::group_by(L2) %>%
    dplyr::mutate(start=abs(start-(max(end)))) %>%
    dplyr::mutate(end=abs(end-max(end))) %>%
    dplyr::mutate(temp=start) %>%
    dplyr::mutate(start=end) %>%
    dplyr::mutate(end=temp) %>%
    dplyr::select(-temp) %>%
    dplyr::ungroup()
  neg_L1L2 <- neg_str %>% dplyr::filter(!(type=="exon") & !(type =="intron") & !(type =="CDS") & !(type=="DELETION"))
  neg_corr <- rbind(neg_L1L2,neg_L3) %>% dplyr::arrange(start)
  
  pos_str <- GFFbyL2 %>% dplyr::filter(strand=="+")
  
  grouped_gff <- rbind(neg_corr,pos_str)
  
  labels <- grouped_gff %>%
    dplyr::filter(type=="gene") %>%
    dplyr::select(tag,ngroup,strand,end)
  
  # CDS_start <- grouped_gff %>% 
  #   dplyr::group_by(L2) %>%
  #   dplyr::filter(type=="CDS") %>% 
  #   dplyr::filter(start==(min(start))) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::filter(start==max(start)) %>%
  #   dplyr::distinct(start,.keep_all = T)
  # 
  # maxCDS <- CDS_start %>% dplyr::select(start,ngroup)
  # 
  # reference_tran <- grouped_gff %>% dplyr::filter(ngroup==maxCDS$ngroup)
  # sub_tran <- grouped_gff %>% 
  #   dplyr::filter(!(ngroup==maxCDS$ngroup)) %>%
  #   dplyr::group_by(L2) %>%
  #   dplyr::mutate(CDSstr=ifelse(type=="CDS",start,1e9)) %>%
  #   dplyr::mutate(hjust_cds=min(CDSstr)) %>%
  #   dplyr::mutate(shift=maxCDS$start-hjust_cds) %>%
  #   dplyr::mutate(start=start+(maxCDS$start-hjust_cds)) %>%
  #   dplyr::mutate(end=end+(maxCDS$start-hjust_cds)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(-CDSstr,-hjust_cds,-shift)
  
  # plottable_df <- rbind(reference_tran,sub_tran)
  plottable_df <- grouped_gff
  
  geneL <- plottable_df %>%
    dplyr::filter(type=="gene") %>%
    dplyr::mutate(len=end-start)
  
  span <- max(geneL$len)
  
  tipE <- plottable_df %>% 
    dplyr::group_by(L2) %>%
    dplyr::filter(type=="exon") %>%
    dplyr::mutate(tip=ifelse(strand=="+",max(start),min(start))) %>%
    dplyr::filter(start==tip) %>%
    dplyr::ungroup()
  
  poly<- tipE %>% 
    group_by(L2) %>%
    dplyr::mutate(xpol=list(c(end/1e3+0.015*(span/1e3),end/1e3+0.015*(span/1e3),end/1e3+0.03*(span/1e3)))) %>%
    dplyr::mutate(ypol=list(c(ngroup+0.05,ngroup-0.05,ngroup))) %>%
    dplyr::mutate(xmin=end/1e3,xmax=end/1e3+0.015*(span/1e3),ymin=ngroup-0.025,ymax=ngroup+0.025)
  
  restE <- plottable_df %>% 
    dplyr::filter(type=="exon")
  
  cdsE <- plottable_df %>% 
    dplyr::filter(type=="CDS")
  
  delE <- plottable_df %>%
    dplyr::filter(type=="DELETION")
  
  ids <- factor(seq(1,nrow(tipE),1))
  t2 <- data.frame(x=unlist(dplyr::pull(poly,xpol)),y=unlist(dplyr::pull(poly,ypol)),z=rep(ids,each=3))
  
  ggplot() + geom_rect(data = restE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25),fill="grey") +
    geom_rect(data = cdsE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.25 , ymax=ngroup+0.25, fill=spp)) +
    geom_rect(data = delE, aes(xmin = start/1e3 ,xmax = end/1e3 ,ymin=ngroup-0.35 , ymax=ngroup-0.40),fill="black") +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=start/1e3,xend=(start+((end-start)/2))/1e3,y=ngroup,yend=ngroup+0.25),color="darkgrey") +
    geom_segment(data = plottable_df %>% dplyr::filter(type=="intron"), aes(x=(start+((end-start)/2))/1e3,xend=end/1e3,y=ngroup+0.25,yend=ngroup),color="darkgrey") +
    #geom_rect(data = poly, aes(xmin = xmin ,xmax = xmax ,ymin=ymin , ymax=ymax),fill="black") +
    #geom_polygon(data = t2,aes(x=x,y=y,group=z),fill="black") +
    annotate('text',x=0,y=labels$ngroup+0.35,label=paste0((labels$tag)),parse=F,fontface='italic',hjust=0,size=8/.pt) +
    #geom_text(data = labels, aes(label=paste0("(",strand,")"),x=(end/1e3)+0.015*(span/1e3),y=ngroup,hjust=0),size=3)+
    theme(axis.title.y  = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          legend.position="none",
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6)) + 
    xlab("Transcript length (kb)") + scale_fill_manual(values=c("#FFA500","lightblue"))
}

getGeneFeature <- function(path,gene,dels) {
  
  gff <- ape::read.gff(file = path) %>% 
    dplyr::mutate(source="AndersenLab")
  
  geneF <- gff %>% 
    dplyr::filter(type =="gene") %>%
    tidyr::separate(attributes, into = c("ID", "post"),sep = ";Name=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("gene:","",ID)) %>%
    dplyr::mutate(ID=gsub("Gene:","",ID)) %>%
    tidyr::separate(post,into=c("pre","aliases"), sep=";Alias=") %>%
    tidyr::separate(aliases,into=c("mAlias","rest"), sep=",",extra = "merge") %>%
    dplyr::mutate(L1=ID)
  
  L1 <- geneF %>% 
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  geneAlias <- geneF %>% 
    dplyr::select(ID,mAlias) 
  
  L2 <- gff %>% 
    dplyr::filter(type=="mRNA") %>%
    tidyr::separate(attributes,into = c("ID","post"),sep=";Parent=") %>%
    dplyr::mutate(ID=gsub("ID=","",ID)) %>%
    dplyr::mutate(ID=gsub("transcript:","",ID)) %>%
    dplyr::mutate(ID=gsub("Transcript:","",ID)) %>%
    tidyr::separate(post,into = c("L1","rest"),sep=";Name=") %>%
    dplyr::mutate(L1=gsub("gene:","",L1)) %>%
    dplyr::mutate(L1=gsub("Gene:","",L1)) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  parentL2 <- L2 %>% 
    dplyr::select(ID,L1)
  
  children <- gff %>% 
    dplyr::filter(!(type=="mRNA") & !(type=="gene")) %>%
    dplyr::filter(!(type=="three_prime_UTR") & !(type=="five_prime_UTR")) %>%
    tidyr::separate(attributes,into=c("L3","L2"),sep=";") %>%
    dplyr::mutate(L2=ifelse(is.na(L2) | grepl("Note=",L2),L3,L2)) %>%
    dplyr::mutate(L3=gsub("ID=","",L3)) %>%
    dplyr::mutate(L3=gsub("CDS:","",L3)) %>%
    dplyr::mutate(L2=gsub("Parent=","",L2)) %>%
    dplyr::mutate(L2=gsub("transcript:","",L2)) %>%
    dplyr::mutate(L2=gsub("Transcript:","",L2)) %>%
    dplyr::mutate(L3=gsub("Parent=","",L3)) %>%
    dplyr::mutate(L3=gsub("transcript:","",L3)) %>%
    dplyr::mutate(L3=gsub("Transcript:","",L3)) 
  
  if(any(grepl(",",children$L2))) {
    children <- children %>%
      tidyr::separate_rows(L2,sep = ",") %>%
      dplyr::mutate(L3=L2) %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) %>%
      dplyr::mutate(L3=paste0(L3,"_",type)) %>%
      dplyr::group_by(L3) %>%
      dplyr::mutate(L3=paste0(L3,row_number()))%>%
      dplyr::ungroup()
  } else {
    children <- children %>%
      dplyr::left_join(parentL2,by = c("L2"="ID")) 
  }
  
  parentL3 <- children %>% 
    dplyr::select(L2,L3)
  
  L3 <- children %>% 
    dplyr::rename(ID=L3) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1)
  
  tabGFF <- BiocGenerics::rbind(L1,L2,L3) %>%
    dplyr::arrange(seqid,start) %>%
    dplyr::filter(L1==gene)
  
  geneMain <- tabGFF %>% 
    dplyr::filter(type=="gene")
  
  subHead <- tabGFF %>% 
    dplyr::filter(type=="mRNA")
  
  geneHead <- geneMain %>%
    dplyr::mutate(L2=paste(subHead$ID,collapse = ",")) %>%
    tidyr::separate_rows(L2,sep = ",")
  
  gFeatures <- tabGFF %>% dplyr::filter(!(type=="gene")) %>%
    dplyr::left_join(parentL3,by=c("ID"="L3")) %>%
    dplyr::mutate(L2=ifelse(type=="mRNA",ID,L2)) 
  
  delFeatures <- dels %>%
    dplyr::left_join(geneAlias,by = 'mAlias') %>%
    dplyr::filter(ID==gene) %>%
    dplyr::rename(L1=ID) %>%
    dplyr::mutate(type="DELETION",strand=unique(gFeatures$strand)) %>%
    dplyr::left_join(parentL2,by="L1") %>%
    dplyr::rename(L2=ID) %>%
    dplyr::mutate(ID=paste(L2,type,sep="_")) %>%
    dplyr::select(seqid,start,end,type,strand,ID,L1,L2)
  
  remFeatures <- rbind(gFeatures,delFeatures)
  
  if (nrow(remFeatures %>% dplyr::filter(type=="intron")) == 0) {
    introns <- remFeatures %>% 
      dplyr::filter(type=="exon") %>% 
      dplyr::group_by(L2) %>%
      dplyr::mutate(istart=end+1,iend=lead(start)-1) %>%
      dplyr::select(-start,-end) %>%
      dplyr::rename(start=istart,end=iend) %>%
      dplyr::filter(!is.na(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(type="intron") %>%
      dplyr::select(seqid,start,end,type,strand,ID,L1,L2) %>%
      dplyr::mutate(ID=gsub("exon","intron",ID))
    
    finalFeatures <- rbind(geneHead,remFeatures,introns) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      dplyr::left_join(geneAlias,by = c("L1"="ID"))%>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start)) 
    
  } else {
    finalFeatures <- rbind(geneHead,remFeatures) %>%
      dplyr::group_by(L2) %>%
      dplyr::arrange(start) %>%
      dplyr::left_join(geneAlias,by = c("L1"="ID")) %>%
      dplyr::mutate(startPos=min(start)) %>%
      dplyr::mutate(end=end-min(start)) %>%
      dplyr::mutate(start=start-min(start))
    
  }
  
  
  features <- dplyr::group_split(finalFeatures)
  
  return(features)
}


path <- "/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS283/N2.WBonly.WS283.PConly.gff3"
dels <- readr::read_tsv("/projects/b1059/projects/Nicolas/collabs/forAmanda/glc-del.tsv")


gene <- "WBGene00000248"
features <- getGeneFeature(path,gene,dels)
gene <- "WBGene00000232"
features2 <- getGeneFeature(path,gene,dels)
gene <-"WBGene00000233" 
features3 <- getGeneFeature(path,gene,dels)
gene <-"WBGene00001591"
features4 <- getGeneFeature(path,gene,dels)


#plotNormGeneFeatures(features[[1]] %>% dplyr::mutate(tag="N2 ben-1") %>% dplyr::mutate(spp="N2"))


sub_gff <- rbind(features[[1]] %>% dplyr::mutate(tag="ben-1") %>% dplyr::mutate(spp="N2"),
                 features2[[4]] %>% dplyr::mutate(tag="avr-14") %>% dplyr::mutate(spp="N2"),
                 features3[[1]] %>% dplyr::mutate(tag="avr-15") %>% dplyr::mutate(spp="N2"),
                 features4[[1]] %>% dplyr::mutate(tag="glc-1") %>% dplyr::mutate(spp="N2"))

p1 <- plotNormGeneFeatures(sub_gff,T)


ggsave(plot = p1,"/projects/b1059/projects/Nicolas/collabs/forAmanda/glc_genes_wDel.tiff",device = 'tiff',dpi = 600,width = 4,height = 2.5,units = 'in')


####### BEN-1 ORTHOGROUP
# gene <- "WBGene00000248"
# features <- getGeneFeature(path,gene)
# gene <- "WBGene00006536"
# features2 <- getGeneFeature(path,gene)
# gene <-"WBGene00006537" 
# features3 <- getGeneFeature(path,gene)
# gene <-"WBGene00003171"
# features4 <- getGeneFeature(path,gene)
# gene <-"WBGene00006538" 
# features5 <- getGeneFeature(path,gene)
# gene <-"WBGene00006539" 
# features6 <- getGeneFeature(path,gene)