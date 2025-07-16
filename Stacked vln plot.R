
#--------------------Stacked vln plot---------------------------------------------------------



markergenes_up <- c("gene1","gene2","gene3","gene4","gene5") 



b <- VlnPlot(combined2, markergenes_up, stack = TRUE, sort = FALSE, flip = TRUE, 
             fill.by = "ident", cols = cluster_colors, adjust = 3) +
  theme_classic()+
  theme(legend.position = "right",
        scale_y_discrete(position = "left"),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, face = "bold"),
        strip.background = element_blank(),
        strip.placement = "left"
  )+  
  ylab("") +  
  ggtitle("Upregulated DEG")

b