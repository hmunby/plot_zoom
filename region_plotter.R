#!/usr/bin/env Rscript

# Plot zoomed plot around GWAS hit and highlight specific SNPs
# manhattan_zoom.R [GEMMA association output file] [PLINK ld output file] [focal snp: chr] [focal snp: pos] [test statistic pval] [output prefix]

#### Import arguments

args <- commandArgs(trailingOnly=TRUE)
data <- args[1]
ld_file <- args[2]
testcol <- as.numeric(args[3])
chrcol <- as.numeric(args[4])
pscol <- as.numeric(args[5])
focal_chrom <- as.numeric(args[6])
positions <- as.numeric(unlist((strsplit(args[7],","))))
focal_snp <- positions[1]
sig_line <- as.numeric(args[8])
transform <- as.logical(as.numeric(args[9]))
head <- as.logical(as.numeric(args[10]))
output <- args[11]
gtf_file <- args[12]
win_start <- as.numeric(args[13])
win_end <- as.numeric(args[14])

# Read file 
res <- read.table(data, header=head)
ld <- read.table(ld_file, header=T)
gtf_ensembl <-  rtracklayer::import(gtf_file)

#### Extract points to plot 
# Restrict to single chromosome & sort by position
res <- res[which(res[,chrcol]==focal_chrom),c(chrcol,pscol,testcol)]
res <- res[order(res[,2]),]

# Check if focal_snp is present in plotting data
if( !(focal_snp %in% res[,2]) ){
  stop("Focal SNP is not present in input data")
} 

# Check if focal_snp is present in LD input 
if( !(focal_snp %in% ld$BP_B) ){
  stop("Focal SNP is not present in VCF")
}

# Match SNPs from two input files, combine and subset to window
res <- res[which(res[,2] %in% ld$BP_B),]
ld <- ld[which(ld$BP_B %in% res[,2]),]
points <- cbind(res, ld$R2)
points <- points[which((points[,2] >= win_start) & (points[,2] <= win_end)),]

#### Specify colours based on LD with focal SNP 
rbPal <- colorRampPalette(c("darkblue","aquamarine","lawngreen","yellow","orange","red"))
points$colour <- rbPal(10)[as.numeric(cut(points$`ld$R2`, breaks = 10))]
fs <- which(points[,2] == focal_snp)
points[fs,"colour"] <- "magenta"

# Convert gtf to data.frame 
gtf_ensembl_df <- as.data.frame(gtf_ensembl)

# Subset to chromosome & region of interest 
gtf_subset_df <- gtf_ensembl_df[which(gtf_ensembl_df$seqnames == focal_chrom),]
gtf_subset_df <- gtf_subset_df[which((gtf_subset_df$start >= win_start & gtf_subset_df$start <= win_end) | (gtf_subset_df$end >= win_start & gtf_subset_df$end <= win_end)),]

# Get list of genes in region
subset_genes <- gtf_subset_df[which(gtf_subset_df$type == "gene"),c("gene_id","start","end")]

# Calculate gene levels
margin <- 1000
subset_genes$level <- c()
subset_genes$t_start <- c()
subset_genes$t_end <- c()
gtf_subset_df$level <- c()
for (i in 1:length(subset_genes$gene_id)){
  # calculate extent of transcripts 
  subset_genes[i,"t_start"] <- min(gtf_subset_df[which((gtf_subset_df$gene_id == subset_genes[i,"gene_id"]) & (gtf_subset_df$type == "transcript")),"start"])
  subset_genes[i,"t_end"] <- max(gtf_subset_df[which((gtf_subset_df$gene_id == subset_genes[i,"gene_id"]) & (gtf_subset_df$type == "transcript")),"end"])
  if(i == 1){
    subset_genes[i,"level"] <- 1
  }else{
    if(subset_genes[i,"t_start"] <= (subset_genes[(i-1),"t_end"]+margin)){
      subset_genes[i,"level"] <- subset_genes[(i-1),"level"] + 1
      if(subset_genes[i,"level"] >= 3){
        last_entry_n_2 <- max(which(subset_genes$level == (subset_genes[i,"level"]-2)))
        if(subset_genes[i,"t_start"] >= subset_genes[last_entry_n_2,"t_end"]+margin){
          subset_genes[i,"level"] <- subset_genes[i,"level"] - 2 
        }
      }
    }else{
      subset_genes[i,"level"] <- 1
    }
  }
  gtf_subset_df[which(gtf_subset_df$gene_id == subset_genes[i,"gene_id"]), "level"] <- subset_genes[i,"level"]
}
levels <- max(gtf_subset_df$level, na.rm = TRUE)

##  Function to plot exon annotation (subplot)
plot_exon <- function(exon_row, b1, t1){
  rect(xleft = exon_row["start"], xright = exon_row["end"], ybottom=b1, ytop=t1, col = "red" , border = NA)
}

## Function to plot genes (subplot)
plot_gene <- function(gene_id, gtf, levels){
  # subset the gtf to only gene of interest
  gtf_mygene <- gtf[which(gtf$gene_id == gene_id),]
  # Assign position based on levels
  gene_level <- gtf_mygene[1,"level"]
  level_width <- 1/levels
  m <- level_width*(levels-gene_level)+0.5*level_width
  t1 <- m + 0.3*level_width
  b1 <- m - 0.3*level_width
  # get label 
  if(!is.na(gtf_mygene[1,"gene_name"])){
    lab <- gtf_mygene[1,"gene_name"]
  } else{
    lab <- gene_id
  }
  if(gtf_mygene[1,"strand"] == "+"){
    lab <- paste0(lab," >")
  }else{
    lab <- paste0("< ",lab)
  }
  # subset exons 
  exons_gtf <- gtf_mygene[which(gtf_mygene$type == "exon"),]
  # plot exons 
  apply(exons_gtf, 1, plot_exon, b1=b1, t1=t1)
  # plot extent of transcripts
  segments(x0 = subset_genes[which(subset_genes$gene_id == gene_id), "t_start"], y0 = m, x1 = subset_genes[which(subset_genes$gene_id == gene_id), "t_end"], y1 = m ,col = "red", lty = 1, lwd=3) 
  # plot 3UTR and 5UTR if present 
  UTR_3 <- which(gtf_mygene$type == "three_prime_utr")
  UTR_5 <- which(gtf_mygene$type == "five_prime_utr")
  if(length(UTR_3) != 0){
    rect(xleft = gtf_mygene[UTR_3,"start"], xright = gtf_mygene[UTR_3,"end"], ybottom=b1, ytop=t1, col = NA , border = "red", lwd = 0.1)
  }
  if(length(UTR_5) != 0){
    rect(xleft = gtf_mygene[UTR_5,"start"], xright = gtf_mygene[UTR_5,"end"], ybottom=b1, ytop=t1, col = NA , border = "red", lwd = 0.1)
  }
  # plot labels 
  midpoint = (gtf_mygene[1,"start"]+gtf_mygene[1,"end"])/2
  text(x=midpoint, y=t1+0.2*level_width, labels=lab, cex = 1.5*level_width, col = "red")
}


## Function to plot SNP density panel
plot_SNPs <- function(SNP_position_vector){
  plot(x=SNP_position_vector, y=rep(1,length(SNP_position_vector)), type="h", xlim=c(win_start,win_end), ylim=c(0,1), axes=FALSE, xlab="", ylab="SNPs", lwd=1.5)
}

## Function to add position labels 
plot_label <- function(position, pos_yval_colour){
  row <- which(pos_yval_colour[,1] == position)
  text(x=pos_yval_colour[row,1], y=(pos_yval_colour[row,2]+0.03*max(pos_yval_colour[,2])), labels=paste0("chr",focal_chrom,":",position), cex = 1.2)
}

## Function to plot test values panel
plot_test_values <- function(pos_yval_colour, log_transform, sig_threshold){
  
  if(log_transform == TRUE){
    # Test if there are any zero
    if(pos_yval_colour[fs,2] == 0){
      stop("Log transforming requested when focal SNP has value zero.")
    }
    if(0 %in% pos_yval_colour[,2]){
      zero_vals <- sum(pos_yval_colour[,2] == 0)
      print(paste0("Warning: Cannot plot log transform for ",zero_vals," SNP(s) as they have test value zero."))
      pos_yval_colour <- pos_yval_colour[-c(which(pos_yval_colour[,2] == 0)),]  # remove any rows with test value 0 from data
    }
    # Log transform yvals
    pos_yval_colour[,2] <- -log(pos_yval_colour[,2], base=10)
    # Log transform significance threshold 
    sig_threshold <- -log(sig_threshold, base=10)
    # Set ylims
    minval <- min(pos_yval_colour[,2])
    maxval <- max(pos_yval_colour[,2])
    if(minval >= 0){
      miny <- 0
    }else{
      miny <- minval - 0.02*abs(minval)
    }
    if(maxval <= 0){
      maxy <- 0
    }else{
      maxy <- maxval + 0.02*maxval
    }
    # Set ylabel 
    labely <- expression("-log"[10] * "(p)")
  }else{
    ## For non-log transformed case 
    # Set ylims based on data spread
    st <- sd(pos_yval_colour[,2])
    miny <- min(pos_yval_colour[,2]) - st
    maxy <- max(pos_yval_colour[,2]) + st
    # Set ylabel 
    if(head == TRUE){
      labely <- colnames(pos_yval_colour[2])
    }else{
      labely <- "Test Value"
    }
  }
  # Remove any rows with NA or non-finite values in test column
  pos_yval_colour <- pos_yval_colour[c(which(is.finite(pos_yval_colour[,2]))),]
  
  # Plotting
  plot(x=pos_yval_colour[,1], y=pos_yval_colour[,2], pch=16, col=pos_yval_colour[,3], xlab=paste0("Position on chromosome ",focal_chrom), ylab=labely, cex.lab=1.2, xlim=c(win_start,win_end), ylim=c(miny,maxy), xaxt='n')
  legend("topleft", title=expression("R"^2),legend=c("focal SNP",1,0.8,0.6,0.4,0.2,0), col=c("magenta", rev(rbPal(6))), pch=20, cex=1.2,  bty="n")
  # Add labels to focal SNP and other specified SNPs
  lapply(positions, plot_label, pos_yval_colour=pos_yval_colour)
  # Add significance threshold to plot if specified 
  if(!is.na(sig_threshold)){
    abline(h=sig_threshold, col = "red", lwd=3) 
  }
}

## Function to plot annotations panel
plot_annotations <- function(gene_id, gtf_subset_df){
  plot(0,type='n', axes = FALSE, ann=FALSE, xlim=c(win_start, win_end), ylim=c(0, 1.15), frame.plot= TRUE)
  axis(3, xpd = TRUE)
  lapply(gene_id, plot_gene, gtf=gtf_subset_df, levels=levels)
}

# Set annotation panel height
annot_h = 0.25+0.75*levels
total_h = 9+annot_h
png_h = 210*total_h

# Plot all panels
png(filename=paste0(output,"_regional_plot.png"), width=2800, height=png_h)
plot.new()
layout(matrix(c(1,0,2,0,3),ncol=1), widths=c(15,15,15,15,15), heights=c(0.8,0.6,6,0.6,annot_h))
annot_h 
par(mar = c(0,5,0,0), oma = c(5, 3, 2, 2), cex=2.5)
plot_SNPs(points[,2])
plot_test_values(points[,c(2,3,5)], transform, sig_line)
plot_annotations(subset_genes[,"gene_id"], gtf_subset_df)
mtext(text=paste0("Position on chromosome ", focal_chrom),side=1,line=1,outer=TRUE,cex=4)
invisible(dev.off())

# Ouput list of SNPs in window and their LD values
colnames(points) <- c("chr","pos","test_value","R2_with_focalSNP","colour")
write.table(points[,c(1,2,3,4)], paste0(output,"_plotted_SNP_values.txt"), quote=FALSE, sep = "\t", col.names=TRUE, row.names = FALSE)
