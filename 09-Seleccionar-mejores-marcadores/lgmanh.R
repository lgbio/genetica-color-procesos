library (ggrepel) # For SNP labels
lgmanhTraits <- function (markers) {
	gwasResults = markers %>% rename ("BP"="POS") 

	# First of all, we need to compute the cumulative position of SNP.
	MAXBP = max (gwasResults$BP)
	don <- gwasResults %>% 
	  # Compute chromosome size
	  group_by(CHR) %>% 
	  summarise(chr_len=max(BP)) %>% 
	  #summarise(chr_len=MAXBP) %>% 
	  # Calculate cumulative position of each chromosome
	  mutate(tot=cumsum(chr_len)-chr_len) %>%
	  select(-chr_len) %>%
	  # Add this info to the initial dataset
	  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
	  # Add a cumulative position of each SNP
	  arrange(CHR, BP) %>%
	  mutate( BPcum=BP+tot)
  	view (don)
	
	# Then we need to prepare the X axis. Indeed we do not want to display 
	# the cumulative position of SNP in bp, but just show the chromosome name instead.
	axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

	# Handle trait names for y axis
	don = don %>% group_by (TRAIT) %>% mutate (gid=cur_group_id()); view (don)

	# Add labels
  	don = don %>% mutate (label="yes")

	# Ready to make the plot using ggplot2:
	ggplot(don, aes(x=BPcum, y=TRAIT)) +
	    # Show all points
	    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.5) +
	    scale_color_manual(values = rep(c("orange", "midnightblue"), 22 )) +
	    # custom X axis:
	    scale_x_continuous (label = axisdf$CHR, breaks= axisdf$center ) +
	    #scale_y_continuous (expand = c(0, 0) ) +     # remove space between plot area and x axis
	    # SNP labels
	    geom_label_repel (data=subset (don, label=="yes"), aes (label=SNP), size=2) +
	    # Custom the theme:
	    theme_bw() + theme(
	      legend.position="none",
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.x = element_blank(),
	      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
	    ) +
        labs (title="Significant SNPs by trait", y="TRAITS", x="CHROMOSOME")
}
