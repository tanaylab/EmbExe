

atlas_proj = function(mat_query,
                      ref_mc_id,
                      ref_mc2d_id,
                      fn,
                      feat_genes,
                      cex.points = 0.5,
                      w = 1000,h = 1000,
                      plot_pdf = F,
                      main_tag = "",
                      cex.main = 1) {
  
  reg = 1e-5
  mc_wt = scdb_mc(ref_mc_id)
  mc2d_wt = scdb_mc2d(ref_mc2d_id)

  feat_genes_f = intersect(feat_genes,intersect(rownames(mat_query),rownames(mc_wt@e_gc)))
  
  message("The following feature genes were excluded:")
  print(setdiff(feat_genes,feat_genes_f))
  
  compute_projection = project_query_umi_matrix_on_metacell_egc_matrix(mat_query = mat_query[feat_genes_f,],
                                                                       mc_egc = mc_wt@e_gc[feat_genes_f,],
                                                                       reg)
  
  best_reference_metacell = compute_projection$best_reference_metacell
  
  
  plot_projected_cells_on_mc2d_background(mc2d_id = ref_mc2d_id,
                                          best_reference_metacell = best_reference_metacell,
                                          filename = fn,
                                          cex.points = cex.points,
                                          cex.main = cex.main)
 

  
  names(best_reference_metacell) = colnames(mat_query)
  
  return(best_reference_metacell)
}

plot_projected_cells_on_mc2d_background = function(mc2d_id,         
                                                   best_reference_metacell,
                                                   filename = NULL,        
                                                   atlas_cells_background_color = NULL,
                                                   cex.points = 0.5,
                                                   w = 1000,
                                                   h = 1000, 
                                                   plot_pdf = F,
                                                   main_tag = "",
                                                   cex.main = 1) {
  
  mc2d_wt = scdb_mc2d(mc2d_id)
  mc_wt = scdb_mc(mc2d_wt@mc_id)
  if (is.null(atlas_cells_background_color)) {
    
    atlas_cells_background_color = rep("gray95",length(mc2d_wt@sc_x))
    names(atlas_cells_background_color) = names(mc2d_wt@sc_x)
    
  }
  
  n = length(best_reference_metacell)
  
  xrange = 0.02 * (max(mc2d_wt@mc_x) - min(mc2d_wt@mc_x))
  yrange = 0.02 * (max(mc2d_wt@mc_y) - min(mc2d_wt@mc_y))
  ref_x = mc2d_wt@mc_x[best_reference_metacell] + rnorm(n, 0, xrange)
  ref_y = mc2d_wt@mc_y[best_reference_metacell] + rnorm(n, 0, yrange)
  xlim = c(min(mc2d_wt@mc_x), max(mc2d_wt@mc_x))
  ylim = c(min(mc2d_wt@mc_y), max(mc2d_wt@mc_y))

  if(!is.null(filename)) {
    if(plot_pdf) {
      pdf(file = filename,width = w,height = h,useDingbats = F)
    } else {
    png(filename = filename,width = w,height = h)
  }
  
  }
 
  # filter cells wit NAs
  f = !is.na(mc2d_wt@sc_x)
  all_cells = names(mc2d_wt@sc_x)[f]
  
  plot(mc2d_wt@sc_x[all_cells],mc2d_wt@sc_y[all_cells], col = atlas_cells_background_color[all_cells],pch = 19,cex = cex.points,xaxt = 'n',yaxt = 'n',xlab = "",ylab = "",axes = F,main = main_tag,
       cex.main = cex.main)
  points(ref_x, ref_y, pch = 19, col = mc_wt@colors[best_reference_metacell], 
         ylim = ylim, xlim = xlim,cex = cex.points)

  if(!is.null(filename)) {
      dev.off()
  }
  
  
}


project_query_umi_matrix_on_metacell_egc_matrix = function(mat_query,mc_egc,reg = 1e-5) {
  
  common_genes = intersect(rownames(mat_query),rownames(mc_egc))
  
  legc = log2(mc_egc[common_genes,] + reg)
  
  query.reference.correlation.matrix = tgstat::tgs_cor(as.matrix(log2(mat_query[common_genes,] + 1)),legc[common_genes,])
  best_ref = apply(query.reference.correlation.matrix,1,which.max)
  
  
  return(list(query.reference.correlation.matrix = query.reference.correlation.matrix,best_reference_metacell = best_ref))
} 


project_query_mc_egc_on_reference_mc_egc = function(egc_query,egc_reference,reg = 1e-5) {
  
  common_genes = intersect(rownames(egc_query),rownames(egc_reference))
  
  legc_reference = log2(egc_reference[common_genes,] + reg)
  legc_query = log2(egc_query[common_genes,] + reg)
  
  query.reference.correlation.matrix = tgstat::tgs_cor(legc_query,legc_reference)
  best_ref = apply(query.reference.correlation.matrix,1,which.max)
  
  return(list(query.reference.correlation.matrix = query.reference.correlation.matrix,best_reference_metacell = best_ref))
} 
