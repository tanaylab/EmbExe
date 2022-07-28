# Plot 2d projection with ggplot
manifold2dproj_hr <- function(mc2d_id = NULL, mc_id = NULL, mat_id = NULL, 
                              plot_dir = NULL, plot_nm = NULL, 
                              plot_mc = F,focus_cells = NULL, plot_par = NULL, plot_title=NULL, 
                              use_mc_egc = F,studio_plot=F,edge_width = 0.2, 
                              sc_point_size = 1, mc_txt_col = NULL,mc_stroke_size = 0.1,
                              mc_point_size = 5, mc_name_txt_size = 2, sc_transparency = 1,
                              pdf_width = 8, pdf_height = 8,plot_gene=NULL, 
                              plot_png = F, png_width = 800,
                              png_height = 800, png_res = 300){
  if(is.null(mc2d_id)){stop("please specify a mc2d ID <missing par: mc2d_id>")}
  if(is.null(mc_id)){stop("please specify a mc ID <missing par: mc_id>")}
  
  
  fn <- sprintf("%s%s",plot_dir,plot_nm)
  mc2d <- scdb_mc2d(mc2d_id)
  mc <- scdb_mc(mc_id)
  
  mc2d_df_sc <- data.frame(sc_x = mc2d@sc_x[names(mc@mc)], 
                           sc_y = mc2d@sc_y[names(mc@mc)], 
                           sc_col = as.character(mc@colors[mc@mc]),
                           stringsAsFactors = F)
  
  mc2d_df_mc <- data.frame(mc_x = mc2d@mc_x, 
                           mc_y = mc2d@mc_y, 
                           mc_col = as.character(mc@colors),
                           stringsAsFactors = F)
  
  mc2d_df_sc <- mc2d_df_sc[!is.na(mc2d_df_sc$sc_x),]
  mc2d_df_mc <- mc2d_df_mc[!is.na(mc2d_df_mc$mc_x),]
  
  edge_start <- mc2d@graph$mc1 # vector coordinates start
  edge_end <- mc2d@graph$mc2 # vector coordinates end
  dx = mc2d@mc_x[edge_start] - mc2d@mc_x[edge_end] # distance over x
  dy = mc2d@mc_y[edge_start] - mc2d@mc_y[edge_end] # distance over y
  f = sqrt(dx^2 + dy^2) > 0 # pitagorian vector size (filtering vectors with size = 0)
  
  mc2d_df_edge <- data.frame(st_x = mc2d@mc_x[edge_start], 
                             st_y = mc2d@mc_y[edge_start],
                             en_x = mc2d@mc_x[edge_end], 
                             en_y = mc2d@mc_y[edge_end])
  
  mc2d_theme <- theme(panel.grid = element_blank(), 
                      panel.border = element_rect(fill=NA),
                      axis.title = element_blank(), 
                      axis.text = element_blank(), 
                      axis.ticks = element_blank(), 
                      axis.line.x = element_blank(), 
                      axis.line.y= element_blank(), 
                      legend.position="None")   
  
  sc <- ggplot(mc2d_df_sc, aes(sc_x,sc_y)) + 
    geom_point(color = mc2d_df_sc$sc_col, size = sc_point_size, alpha = sc_transparency) + 
    scale_color_identity() + mc2d_theme
  
  if(!is.null(plot_gene)){
    gene <- rownames(mc@mc_fp)[which(rownames(mc@mc_fp)==plot_gene)]
    if(! gene %in% rownames(mc@mc_fp)){}
    lfp <- log2(mc@mc_fp[gene, ])
    max_lfp <- max(lfp)
    min_lfp <- min(lfp)
    lfp_n =lfp
    lab = "Fold Change"
    lfp_n = pmin(pmax(lfp, min_lfp), max_lfp) - min_lfp
    if(use_mc_egc){
      lfp <- log2(mc@e_gc[gene, ] + 1e-5)
      max_lfp <- max(lfp)
      min_lfp <- min(lfp)
      lfp_n = lfp
      lab = "Log2(Absolute expression)"}
    
    mc2d_df_mc <- data.frame(mc_x = mc2d@mc_x, 
                             mc_y = mc2d@mc_y, 
                             gene_expression = lfp_n,
                             stringsAsFactors = F)
    
    mc2d_df_mc <- mc2d_df_mc[!is.na(mc2d_df_mc$mc_x),]
    
    mc2d_gene_theme <- theme(panel.grid = element_blank(), 
                        panel.border = element_rect(fill=NA),
                        axis.title = element_blank(), 
                        axis.text = element_blank(), 
                        axis.ticks = element_blank(), 
                        axis.line.x = element_blank(), 
                        axis.line.y= element_blank())
    
    sc1 <- ggplot(data = mc2d_df_mc, aes(mc_x,mc_y,fill=gene_expression)) + 
      mc2d_gene_theme

    mc2d_df_mc$ct <- annotation_col
    
    sum <- mc2d_df_mc %>% 
      group_by(ct) %>% 
      summarise(mean = mean(gene_expression),sd=sd(gene_expression))
    
    ct_high <- unlist(sum[which.max(sum$mean),"ct"])[[1]]
    low_level <- colorRampPalette(RColorBrewer::brewer.pal(9,"Greys"))(1000)[1]
    high_level <- tail(colorRampPalette(RColorBrewer::brewer.pal(9,"Greys"))(1000))[1]
    
    sc <- sc1 + geom_point(data = mc2d_df_sc, aes(sc_x,sc_y),
                 shape=19, fill = "white", color = mc2d_df_sc$sc_col,size = sc_point_size) + 
      geom_point(size = mc_point_size, pch = 21, stroke=mc_stroke_size, color = "black") + 
      scale_fill_gradient(name = lab,low = low_level,high = high_level) +
      theme(plot.title = element_text(face = "bold",hjust = 0.5), 
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "bottom")  +
      ggtitle(label = sprintf("%s",gene),
              subtitle = sprintf("Highest Expression: %s",ct_high))
  }
  
  if(!is.null(focus_cells)){
    mc2d_df_sc_focus <- data.frame(sc_x = mc2d@sc_x[focus_cells], 
                                   sc_y = mc2d@sc_y[focus_cells], 
                                   sc_col = as.character(mc@colors[mc@mc[focus_cells]]),
                                   stringsAsFactors = F)
    
    mc2d_df_sc_grey <- data.frame(sc_x = mc2d@sc_x[names(mc@mc)], 
                                  sc_y = mc2d@sc_y[names(mc@mc)], 
                                  sc_col = as.character(mc@colors[mc@mc]),
                                  stringsAsFactors = F)
    
    mc2d_df_sc_focus <- mc2d_df_sc_focus[!is.na(mc2d_df_sc_focus$sc_x),]
    mc2d_df_sc_grey <- mc2d_df_sc_grey[!is.na(mc2d_df_sc_grey$sc_x),]
    
    sc_grey <- ggplot(mc2d_df_sc_grey, aes(sc_x,sc_y)) + 
      geom_point(color = "#E5E5E5", size = sc_point_size) + mc2d_theme
    
    sc <- sc_grey + geom_point(data = mc2d_df_sc_focus, aes(sc_x,sc_y),
                               color = mc2d_df_sc_focus$sc_col, size = sc_point_size)
  }
  
  if(!is.null(plot_par)){
    mat <- scdb_mat(mat_id)
    if(is.null(mat_id)){stop("please specify a mat ID <missing par: mat_id>")}
    cells_split <- split(names(mc@mc), mat@cell_metadata[names(mc@mc),plot_par])
    plot_par_levels <- names(cells_split)
    for(level in plot_par_levels){
      cells_sub <- cells_split[[level]]
      cells_sub <- cells_sub[cells_sub %in% names(mc2d@sc_x)]
      fn_lev <- sprintf("%s_%s.pdf",fn,level)
      mc2d_df_sc_focus <- data.frame(sc_x = mc2d@sc_x[cells_sub], 
                                     sc_y = mc2d@sc_y[cells_sub], 
                                     sc_col = as.character(mc@colors[mc@mc[cells_sub]]),
                                     stringsAsFactors = F)
      
      mc2d_df_sc_grey <- data.frame(sc_x = mc2d@sc_x[names(mc@mc)], 
                                    sc_y = mc2d@sc_y[names(mc@mc)], 
                                    sc_col = as.character(mc@colors[mc@mc]),
                                    stringsAsFactors = F)
      
      mc2d_df_sc_focus <- mc2d_df_sc_focus[!is.na(mc2d_df_sc_focus$sc_x),]
      mc2d_df_sc_grey <- mc2d_df_sc_grey[!is.na(mc2d_df_sc_grey$sc_x),]
      
      sc_grey <- ggplot(mc2d_df_sc_grey, aes(sc_x,sc_y)) + 
        geom_point(color = "#E5E5E5", size = sc_point_size) + mc2d_theme
      
      sc <- sc_grey + 
        geom_point(data = mc2d_df_sc_focus, aes(sc_x,sc_y), shape = 21, color = "black",stroke = 0.2,
                                 bg = mc2d_df_sc_focus$sc_col, size = sc_point_size) + 
        theme(plot.title = element_text(face = "bold",hjust = 0.95, vjust = -8),
              axis.title = element_blank())  +
        ggtitle(sprintf("%s %s",plot_title,level))
      if(studio_plot){
        print(sc)
      } else{pdf(file = fn_lev, width = pdf_width, height = pdf_height, useDingbats = F)
        print(sc)
        dev.off()}
      
    }
    
  }
  
  if(plot_mc){
    mc2d_df_mc$sub_txt <- "black"
    mc2d_df_mc$sub_txt[mc2d_df_mc$mc_col == "#67000d"] <- "white"
    if(!is.null(mc_txt_col))
    {mc2d_df_mc$sub_txt <- mc_txt_col}
    sc <- sc + geom_segment(data = mc2d_df_edge, aes(x = st_x, y = st_y, xend = en_x, yend = en_y),
                      size = edge_width) +
      geom_point(data = mc2d_df_mc, aes(mc_x,mc_y),
                 fill = mc2d_df_mc$mc_col, size = mc_point_size, pch = 21, color = "black") +
      geom_text(data = mc2d_df_mc, aes(mc_x,mc_y,label = 1:nrow(mc2d_df_mc),fontface = "bold", 
                                       col = sub_txt),size = mc_name_txt_size)
  }
  
  if(plot_png){
    fn <- paste(fn,".png",sep="")
    png(filename = fn, width = png_width, height = png_height, res = png_res)
    print(sc)
    
  }
  
  if(studio_plot){
    print(sc)
    stop_quietly <- function() {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    }
    message("Studio plot completed")
    stop_quietly()
  }
  
  fn <- paste(fn,".pdf",sep="")
  pdf(file = fn, width = pdf_width, height = pdf_height, useDingbats = F)
  print(sc)
  
  stop_quietly <- function() {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  
  message(sprintf("Completed"))
  dev.off()
  stop_quietly()
}

