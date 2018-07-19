# Here you can find the instructions for plotting the mapped mitochonfrial CR sequences onto phylogeny.
# Before you start please calculate the intervals for genomic regions you are going to map.
# Use MITOS web-server (http://mitos.bioinf.uni-leipzig.de/index.py) for either complete genome annotation or the sequences which flank regions adjacent with CR. 
# For CR annotation use... previously published annotation as nobody knows yet what is a biological scence of those up to 6 conserved elements within it.

# setting up your working directory
setwd(choose.dir())

# you will need the following packages
# but you don't really need xlsx if you know how to deal with tables in R
install.packages("xlsx","genoPlotR")

# to load the packages
library(genoPlotR)
library(xlsx)

# you will need to prepare the tables like that (https://github.com/Sturcoal/Mitochondrial-Control-Region/blob/master/AD.csv) for every tip on your phylogeny
# please read the genoPlotR package manual to find out more

# read the tables
AD <- read.xlsx("AD.xls",1)
AD_seg <- dna_seg(AD)
AL <- read.xlsx("AL.xls",1)
AL_seg <- dna_seg(AL)
AM <- read.xlsx("AM.xls",1)
AM_seg <- dna_seg(AM)
HD <- read.xlsx("HD.xls",1)
HD_seg <- dna_seg(HD)
PT <- read.xlsx("PT.xls",1)
PT_seg <- dna_seg(PT)
PC <- read.xlsx("PC.xls",1)
PC_seg <- dna_seg(PC)
PF <- read.xlsx("PF.xls",1)
PF_seg <- dna_seg(PF)
SG <- read.xlsx("SG.xls",1)
SG_seg <- dna_seg(SG)
CJ <- read.xlsx("CJ.xls",1)
CJ_seg <- dna_seg(CJ)
LF <- read.xlsx("LF.xls",1)
LF_seg <- dna_seg(LF)

# prepare the list with data on each annotated sequence 
dna_segs <- list(PT_seg,AM_seg,AL_seg,AD_seg,PF_seg,PC_seg,LF_seg,SG_seg,CJ_seg,HD_seg)
dna_segs[[1]]$gene_type <- "exons"
dna_segs[[2]]$gene_type <- "exons" 
dna_segs[[3]]$gene_type <- "exons" 
dna_segs[[4]]$gene_type <- "exons"
dna_segs[[5]]$gene_type <- "exons"
dna_segs[[6]]$gene_type <- "exons"
dna_segs[[7]]$gene_type <- "exons"
dna_segs[[8]]$gene_type <- "exons"
dna_segs[[9]]$gene_type <- "exons"
dna_segs[[10]]$gene_type <- "exons"

# a creepy way to extract and attach the names to "DNA_seqs" vector
names <- c("PT","AM","AL","AD","PF","PC","LF","SG","CJ","HD")
names(dna_segs) <- names

# phylogenetic tree must be prepared in newick fromat and transformes into object of class phylog using newick2phylog function
Zo_k80_NJ_tree <- newick2phylog("(((((AL:0.00591133,AM:0.00763547)0.8120:0.00346367,AD:0.01008313)1.0000:0.03586823,LT:0.07466133)0.9460:0.01878079,((PC:0.02001232,PF:0.01939655)1.0000:0.03325123,SG:0.06527094)0.8700:0.01062192)0.0000:0.03863916,HD:0.13208128);")
Z_tree <- newick2phylog("(HD:0.21088325,(PT:0.08044520,(((LF:0.09349587,(SG:0.04955483,CJ:0.08538506)82:0.01940017)72:0.01703654,(PF:0.02117793,PC:0.02368433)98:0.04808177)85:0.03694716,((AM:0.00704376,AL:0.00740218)47:0.00523579,AD:0.00874105)100:0.04161217)52:0.01706786)0:0.03504645);")
Z_tree <- newick2phylog("(((PT:0.10403242559705062,((AM:0.007485885607887377,AL:0.007626243016379486):0.005781395800363687,AD:0.012454195787532651):0.04481099942721428):0.016661649974292703,((PF:0.023584205815936343,PC:0.02499369392870654):0.05406792960414841,(LF:0.11106602374788438,(SG:0.05295740130832083,CJ:0.09954401257653378):0.02308741267671227):0.017559401324998725):0.038281143573640875):0.1790513986353381,HD:0.3581027972706762);")

# one more format for plotting options (see below)
Z_tree <- as.hclust(Z_tree)

## creepy comparisons 
df_7 <- data.frame(start1 = c(SG_seg$start[1:13]), end1 = c(SG_seg$end[1:13]),start2 = HD_seg$start, end2 = HD_seg$end)
comp7 <- comparison(df_7)

df_4 <- data.frame(start1 = PC_seg$start, end1 = PC_seg$end,start2 = PF_seg$start, end2 = PF_seg$end)
comp4 <- comparison(df_4)

df_5 <- data.frame(start1 = PF_seg$start, end1 = PF_seg$end,start2 = c(SG_seg$start[2:13]), end2 = c(SG_seg$end[2:13]))
comp5 <- comparison(df_5)

df_2 <- data.frame(start1 = AL_seg$start, end1 = AL_seg$end,start2 = AM_seg$start, end2 = AM_seg$end)
comp2 <- comparison(df_2)

df_3 <- data.frame(start1 = AM_seg$start, end1 = AM_seg$end,start2 = AD_seg$start, end2 = AD_seg$end)
comp3 <- comparison(df_3)

df_1 <- data.frame(start1 = AD_seg$start, end1 = AD_seg$end,start2 = c(LT_seg$start[2:13]), end2 = c(SG_seg$end[2:13]))
comp1 <- comparison(df_1)

df_4 <- data.frame(start1 = c(PT_seg$start[2:13]), end1 = c(PT_seg$end[2:13]),start2 = PC_seg$start, end2 = PC_seg$end)
comp4 <- comparison(df_4)

# list of comparisons
comparisons <- list(comp1,comp2,comp3,comp4,comp5,comp6,comp7)

# I was needed to prepare the colors good-looking mapping and plotting
colfunc <- colorRampPalette(c("black","white"))
co <- colfunc(10)
co

color <- c(co[5],co[5],co[5],co[3:8],co[c(4,5,7)])
color


# attach the colors to comparisons
comparisons[[1]]$col <- co[8]
comparisons[[2]]$col <- co[8]
comparisons[[3]]$col <- co[8]
comparisons[[4]]$col <- co[8]
comparisons[[5]]$col <- co[8]
comparisons[[6]]$col <- co[8]
comparisons[[7]]$col <- co[8]



# plot and have a look
plot_gene_map(dna_segs,tree=Zo_k80_NJ_tree,comparisons=comparisons)

plot_gene_map(dna_segs,tree=Z_tree)

