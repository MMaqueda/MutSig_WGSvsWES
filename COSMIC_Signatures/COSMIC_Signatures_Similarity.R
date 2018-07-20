library(RColorBrewer)

cos_sim_samples_signatures = cos_sim_matrix(cancer_signatures, cancer_signatures)
#plot_cosine_heatmap(cos_sim_samples_signatures,cluster_rows = FALSE)
#col_order = cosmic_order


require(corrplot)
M <- as.data.frame(cos_sim_samples_signatures)
col<- colorRampPalette(c("white", "red"))(8)

pdf("CosSim_COSMIC.pdf",width=12,height=12)
corrplot(as.matrix(M), method="circle", is.corr=FALSE, addCoef.col = "black",
         tl.col="black", tl.srt=45, tl.cex = 0.9,type="upper", col = brewer.pal(n = 10, name = "RdBu")) #Text label color and rotation,)
dev.off()

#,addCoef.col = "black",number.cex = 0.4