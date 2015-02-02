# set working folder to this directory
# setwd("/data/m.gleditzsch/workspace_git/comparative-sequence-logo/")
source("./diffSeqLogo.R");

motif_folder = "CTCF_motifs";
#motif_names = c("GM12878","HeLa-S3","HepG2","HUVEC","K562","MCF7","NHEK","ProgFib");
motif_names = c("HepG2","MCF7","HUVEC","ProgFib");
margin = 0.03;
ratio=16/10;

dim = length(motif_names);

CTCF = list();

#pdf(paste("CTCF_motifs.pdf",sep=""), width=6, height=4.5);
for (name in motif_names) {
    file = paste("./",motif_folder,"/",name,".txt",sep="");
    CTCF[[name]] = as.matrix(read.delim(file,header=F));

    seqLogo(CTCF[[name]])
    
}
#dev.off();


#pdf(paste("CTCF_motifs_comparison_",type,".pdf",sep=""), width=dim*2, height=dim*2/ratio);
    diffLogoTable(CTCF,margin=margin);
#dev.off();
