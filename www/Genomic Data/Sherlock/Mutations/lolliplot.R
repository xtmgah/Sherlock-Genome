set_wd()
libztw()



# load and prepare data---------------------------------------------------------------
load('../RDS/sherlock_data_all.RData')
load('sherlock_maf.RData')

#tmp <- read_csv('../RDS/oncoplot_colors.csv')
tmp <- read_csv('oncoplot_colors.csv')
landscape_colors <- tmp$Color
names(landscape_colors) <- tmp$Name

gff0 = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
gff0 = readRDS(file = gff0)
gff0 <- as_tibble(gff0)


## input paramters
gene <- 'TP53'
group <- NULL # N_A, NU,S_U
minN <- 5



if(is.null(group)){
  samplelist <- sherlock_overall$Tumor_Barcode
}else{
  samplelist <- sherlock_overall %>% filter(SP_Group %in% group)
}

tdata1 <- sherlock_maf %>%
  filter(Hugo_Symbol %in% gene, Tumor_Barcode %in% samplelist,Variant_Classification != "Splice_Site") %>% 
  mutate(AAChange.refGene = if_else(is.na(AAChange.refGene),'Splice_Site',AAChange.refGene)) %>% 
  select(Tumor_Barcode,Variant_Classification,AAChange.refGene) %>% 
  separate_rows(AAChange.refGene,sep = ',') %>% 
  separate(AAChange.refGene,into = c('Gene','TS','Exon','cDNA','AAChange'),sep = ':') %>% 
  mutate(AAPosition = parse_number(str_remove(AAChange,'p.'))) 


tdata2 <- sherlock_maf %>%
  filter(Hugo_Symbol %in% gene, Tumor_Barcode %in% samplelist,Variant_Classification == "Splice_Site") %>% 
  select(Tumor_Barcode,Variant_Classification,GeneDetail.refGene) %>% 
  separate_rows(GeneDetail.refGene,sep = ';') %>% 
  separate(GeneDetail.refGene,into = c('TS','Exon','cDNA'),sep = ':') %>% 
  mutate(Gene = gene) %>% 
  mutate(AAPosition = ceiling(parse_number(str_remove(cDNA,'^c.'))/3)) %>% 
  mutate(AAChange = paste0('X',AAPosition,'_splice'))

tdata0 <- bind_rows(tdata1,tdata2)

tslist <- tdata0 %>% filter(!is.na(TS),TS %in% gff0$refseq.ID) %>% pull(TS) %>% unique()

tslist_select <- tslist[1]
#tslist_select <- "NM_201283"

tdata <- tdata0 %>% filter(TS==tslist_select | is.na(TS)) %>% unique() 

## label 

genefreq <- percent_format(accuracy = 0.1)(length(unique(tdata$Tumor_Barcode))/length(samplelist))

genetitle <- paste0(gene,' mutated in ',genefreq, ' tumors (N=',length(samplelist),')')


# Check the domain --------------------------------------------------------
gff <- gff0 %>% filter(HGNC==gene,refseq.ID==tslist_select) %>% mutate(AAPosition=(End+Start)/2) %>% mutate(n=-0.4)
aasize <- gff$aa.length[1]
dcolors <- c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
gff$color <- dcolors[1:dim(gff)[1]]

tdata_count <- tdata %>% count(AAPosition,AAChange,Variant_Classification)
ymax <- max(tdata_count$n)*1.05


tdata_count2 <- tdata_count %>% filter(n>5)

tdata_count %>% 
  ggplot(aes(x=AAPosition,y=n))+
  geom_segment(aes(x=AAPosition,xend=AAPosition,y=-0.03*ymax,yend=n),size=0.2)+
  geom_point(aes(fill=Variant_Classification),pch=21,size=4)+
  ggrepel::geom_text_repel(
    data = tdata_count2,
    aes(label=AAChange,color=Variant_Classification),
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.2,
    direction    = "x",
    angle        = 90,
    vjust        = 0.5,
    hjust = -1,
    segment.size = 0.1,
    max.iter = 1e4, max.time = 1
  ) +
  geom_rect(aes(xmin=0,xmax=aasize,ymin=-0.06*ymax,ymax=-0.03*ymax),fill="#cccccc",show.legend=FALSE)+
  geom_rect(data=gff,aes(xmin=Start,xmax=End,ymin=-0.07*ymax,ymax=-0.02*ymax),fill=gff$color,show.legend=FALSE)+
  geom_text(data=gff,aes(x=AAPosition,y=-0.045*ymax,label=Label),family = 'Roboto Condensed',fontface = 'plain',size=3.5)+
  scale_fill_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification)+
  scale_color_manual(values =landscape_colors,breaks = tdata_count$Variant_Classification,guide='none')+
  scale_x_continuous(breaks = pretty_breaks(n = 8),limits = c(0,aasize),expand = expansion(mult = c(0.01,0.01)))+
  scale_y_continuous(breaks = pretty_breaks(),limits = c(-0.08*ymax,ymax),expand = expansion(mult = c(0,0)))+
  guides(color='none',fill=guide_legend(nrow=1))+
  theme_ipsum_rc(base_size = 12,axis_title_just = 'm',axis_title_size = 14,grid=FALSE,axis = T,axis_col = 'black',ticks = T)+
  theme(legend.position = 'bottom',legend.text=element_text(size=12))+
  labs(x='',y=paste0('# ',gene,' Mutations'),fill="Mutation Type",subtitle = genetitle)

ggsave(filename = 'tmp.pdf',width = 18,height = 6,device = cairo_pdf)


tdata
