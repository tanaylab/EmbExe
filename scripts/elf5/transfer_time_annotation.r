library("Matrix")
library("qlcMatrix")


elf5_timing = function(mat_nm,excluded_colors = NULL) {
  
  cgraph_id =  mat_nm
  mat = scdb_mat(mat_nm)
  mc_wt = scdb_mc("embexe_recolored")
  
  # Exclude extraembryonic ectoderm and parietal endoderm
  
  load(sprintf("data/%s/color_annotation/cmp_annot.Rda",mat_nm))
  
  f = mat@cell_metadata[colnames(mat@mat),"cell_type"] %in% c("KO","control","unclear","Elf5")
  query_cells = colnames(mat@mat)[f]
  query_cells = intersect(query_cells,names(cmp_annot$query_cls_col))
  
  if(!is.null(excluded_colors)) {
    excluded_cls = query_cells[cmp_annot$query_cls_col[query_cells] %in% excluded_colors]
    query_cells = setdiff(query_cells,excluded_cls)
  } else {
    message("Did not exclude any cell types for timing the embryos.")
  }
  
  n_cells_per_embryo = table(mat@cell_metadata[query_cells,"embryo"])
  
  query_embryos = names(n_cells_per_embryo)[n_cells_per_embryo > 10]
  
  
  wt_cls = intersect(names(mc_wt@mc)[ !(mc_wt@colors[mc_wt@mc] %in% excluded_colors) ],colnames(mat@mat))
  wt_cls = wt_cls[!is.na(mat@cell_metadata[wt_cls,"transcriptional_rank"])]
  atlas_time = mat@cell_metadata[wt_cls,"transcriptional_rank"]
  names(atlas_time) = wt_cls
  
  data_dir = sprintf("data/%s",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  data_dir = sprintf("data/%s/time_match",mat_nm)
  if(!dir.exists(data_dir)) {
    dir.create(data_dir)
  }
  
  
  # first get atlas time distribution
  
  atlas_time_dist = get_atlas_time_dist(atlas_time = atlas_time,graph_id = cgraph_id)
  
  # first timing using control cells only
  query_cls = query_cells[( mat@cell_metadata[query_cells,"cell_type"] %in%  c("KO","unclear","control","Elf5") ) & ( mat@cell_metadata[query_cells,"embryo"] %in%  query_embryos )]
  query_cls_md = mat@cell_metadata[query_cls,"embryo"]
  names(query_cls_md) = query_cls
  
  query_time_dist= get_query_time_dist(query_cls_md = query_cls_md,atlas_time = atlas_time,graph_id = cgraph_id)
  
  time_dist_query = list(atlas_time_dist = atlas_time_dist$atlas_time_dist,
                        query = query_time_dist$query_time_dist)
  
  save(time_dist_query,file = sprintf("%s/time_dist_query.Rda",data_dir))
  
  chim_emb_summary = as.data.frame.matrix(table(mat@cell_metadata[query_cells,"embryo"],mat@cell_metadata[query_cells,"cell_type"]))
  chim_emb_summary$embryo = rownames(chim_emb_summary)
  chim_emb_summary$best_rank_query = NA
  chim_emb_summary[rownames(time_dist_query$query),"best_rank_query"] = time_dist_best_match(atlas_time_dist = time_dist_query$atlas_time_dist,
                                                                                       query_time_dist = time_dist_query$query)
  
  write.table(chim_emb_summary,file = sprintf("%s/time_match_summary.txt",data_dir),sep ="\t",row.names = F)

}

get_query_time_dist = function(query_cls_md,atlas_time,graph_id) {
  
  cgraph = scdb_cgraph(graph_id)
  
  query_cls = intersect(names(query_cls_md),cgraph@nodes)
  atlas_cls = intersect(names(atlas_time),cgraph@nodes)
  atlas_time = atlas_time[atlas_cls]
  query_cls_md = query_cls_md[query_cls]
  
  
  cell_names = c(1:length(cgraph@nodes))
  names(cell_names) = cgraph@nodes
  
  graph = cgraph@edges
  graph$mc1 = as.factor(graph$mc1)
  graph$mc2 = as.factor(graph$mc2)
  levels(graph$mc1) = cell_names[levels(graph$mc1)]
  levels(graph$mc2) = cell_names[levels(graph$mc2)]
  
  # adjacency matrix 
  knn_mat = sparseMatrix(as.numeric(graph$mc1),as.numeric(graph$mc2),x = graph$w)
  colnames(knn_mat) = cgraph@nodes
  rownames(knn_mat) = cgraph@nodes
  
  knn_mat_f = knn_mat[query_cls,atlas_cls]
  
  a = rowMax(X = knn_mat_f,which = T)
  time_match_ind = summary(a$which)
  
  query_time_match = atlas_time[time_match_ind$j]
  names(query_time_match) = query_cls[time_match_ind$i]
  
  query_time_dist = table(query_cls_md[time_match_ind$i],atlas_time[time_match_ind$j])
  tmp = matrix(0,nrow = nrow(query_time_dist),ncol = length(unique(atlas_time)))
  rownames(tmp) = rownames(query_time_dist)
  colnames(tmp) = sort(unique(atlas_time))
  tmp[rownames(query_time_dist),colnames(query_time_dist)] = query_time_dist
  query_time_dist = tmp
  
  
  return(list(query_time_dist = query_time_dist,query_time_match = query_time_match))
}

get_atlas_time_dist = function(query_cls_md,atlas_time,graph_id) {
  
  
  # query_cls_md is a named vector with query_cls as names and the embryo annotation or sth similar as value
  # atlas_time is a named vector giving the atlas_time for each atlas cell
  
  cgraph = scdb_cgraph(graph_id)
  
  atlas_cls = intersect(names(atlas_time),cgraph@nodes)
  atlas_time = atlas_time[atlas_cls]
  
  
  cell_names = c(1:length(cgraph@nodes))
  names(cell_names) = cgraph@nodes
  
  graph = cgraph@edges
  graph$mc1 = as.factor(graph$mc1)
  graph$mc2 = as.factor(graph$mc2)
  levels(graph$mc1) = cell_names[levels(graph$mc1)]
  levels(graph$mc2) = cell_names[levels(graph$mc2)]
  
  knn_mat = sparseMatrix(as.numeric(graph$mc1),as.numeric(graph$mc2),x = graph$w)
  colnames(knn_mat) = cgraph@nodes
  rownames(knn_mat) = cgraph@nodes
  
  knn_mat_f = knn_mat[atlas_cls,atlas_cls]
  for(i in unique(atlas_time)) {
    f = atlas_time == i
    knn_mat_f[f,f] = 0
  }
  
  a = rowMax(X = knn_mat_f,which = T)
  time_match_ind = summary(a$which)
  
  atlas_time_match = atlas_time[time_match_ind$j]
  names(atlas_time_match) = atlas_cls[time_match_ind$i]
  
  q_time_match_tmp = table(atlas_time[time_match_ind$i], atlas_time_match)
  q_time_match = matrix(0,nrow = length(unique(atlas_time)),ncol = length(unique(atlas_time)))
  rownames(q_time_match) = sort(unique(atlas_time))
  colnames(q_time_match) = sort(unique(atlas_time))
  q_time_match[rownames(q_time_match_tmp),colnames(q_time_match_tmp)] = q_time_match_tmp
  
  return(list(atlas_time_dist = q_time_match,atlas_time_match = atlas_time_match))
}

get_query_and_atlas_time_dist = function(query_cls_md,atlas_time,graph_id) {
  
  # query_cls_md is a named vector with query_cls as names and the embryo annotation or sth similar as value
  # atlas_time is a named vector giving the atlas_time per for each atlas cell
  
  cgraph = scdb_cgraph(graph_id)
  
  query_cls = intersect(names(query_cls_md),cgraph@nodes)
  atlas_cls = intersect(names(atlas_time),cgraph@nodes)
  atlas_time = atlas_time[atlas_cls]
  query_cls_md = query_cls_md[query_cls]
  
  
  
  cell_names = c(1:length(cgraph@nodes))
  names(cell_names) = cgraph@nodes
  
  graph = cgraph@edges
  graph$mc1 = as.factor(graph$mc1)
  graph$mc2 = as.factor(graph$mc2)
  levels(graph$mc1) = cell_names[levels(graph$mc1)]
  levels(graph$mc2) = cell_names[levels(graph$mc2)]
  
  knn_mat = sparseMatrix(as.numeric(graph$mc1),as.numeric(graph$mc2),x = graph$w)
  colnames(knn_mat) = cgraph@nodes
  rownames(knn_mat) = cgraph@nodes
  
  knn_mat_f = knn_mat[atlas_cls,atlas_cls]
  for(i in unique(atlas_time)) {
    f = atlas_time == i
    knn_mat_f[f,f] = 0
  }
  
  time_match = apply(knn_mat_f,1,which.max)
  time_match = atlas_time[time_match]
  
  q_time_match_tmp = table(atlas_time, time_match)
  q_time_match = matrix(0,nrow = length(unique(atlas_time)),ncol = length(unique(atlas_time)))
  rownames(q_time_match) = sort(unique(atlas_time))
  colnames(q_time_match) = sort(unique(atlas_time))
  q_time_match[rownames(q_time_match_tmp),colnames(q_time_match_tmp)] = q_time_match_tmp
  
  
  query_time_match = apply(knn_mat[query_cls,atlas_cls],1,which.max)
  query_time_dist = table(query_cls_md,atlas_time[query_time_match])
  tmp = matrix(0,nrow = nrow(query_time_dist),ncol = length(unique(atlas_time)))
  rownames(tmp) = rownames(query_time_dist)
  colnames(tmp) = sort(unique(atlas_time))
  tmp[rownames(query_time_dist),colnames(query_time_dist)] = query_time_dist
  query_time_dist = tmp
  
  
  return(list(atlas_time_dist = q_time_match,query_time_dist = query_time_dist))
}


time_dist_best_match = function(query_time_dist,atlas_time_dist) {
  
  rank_to_time = read.table(file = "data/embexe.transcriptional_rank_developmental_time.tsv",stringsAsFactors = F,h = T,sep = "\t")
  dev_time = rank_to_time$developmental_time
  
  query_ref_cor = cor(t(as.matrix(as.data.frame.matrix(query_time_dist))),
                          t(as.matrix(as.data.frame.matrix(atlas_time_dist))))
  
  
  wt_dev_time = rank_to_time$developmental_time[as.numeric(colnames(query_ref_cor))] 

  best_fit = apply(query_ref_cor,1,function(emb_cor) {
    smooth_spline = smooth.spline(x = wt_dev_time,y = emb_cor,spar = 0.8)
    y = smooth_spline$y
    return(which.max(y))
  })
  
  best_fit = as.numeric(colnames(query_ref_cor))[best_fit]
  return(best_fit)
}



chimera_plot_cor_time_dist_per_emb = function(mat_nm) {
  
  n_cls_min = 19
  rank_to_time = read.table(file = "data/embexe.transcriptional_rank_developmental_time.tsv",stringsAsFactors = F,h = T,sep = "\t")
  dev_time = rank_to_time$developmental_time
  data_dir = sprintf("data/%s/time_match",mat_nm)
  load(sprintf("%s/time_dist_host_ko.Rda",data_dir))
  
  host_dist_all = time_dist_host_ko$host
  ko_dist_all = time_dist_host_ko$ko
  atlas_time_dist = time_dist_host_ko$atlas_time_dist
  
  
  
  host_ref_cor = tgs_cor(t(as.matrix(as.data.frame.matrix(host_dist_all))),
                         t(as.matrix(as.data.frame.matrix(atlas_time_dist))))
  ko_ref_cor = tgs_cor(t(as.matrix(as.data.frame.matrix(ko_dist_all))),
                       t(as.matrix(as.data.frame.matrix(atlas_time_dist))))
  
  chim_time_summary = read.table(file = sprintf("%s/time_match_summary.txt",data_dir),sep ="\t",h = T,stringsAsFactors = F)
  
  f = ( chim_time_summary$host > n_cls_min )
  
  chim_embryos = chim_time_summary$embryo[f]
  
  if(!dir.exists(sprintf("figs/%s",mat_nm))) {
    dir.create(sprintf("figs/%s",mat_nm))
  }
  
  fig_dir = sprintf("figs/%s/time_match",mat_nm)
  if(!dir.exists(fig_dir)) {
    dir.create(fig_dir)
  }
  
  for (chimera in chim_embryos) {
    
    smooth_spline_host = smooth.spline(x = dev_time,y = host_ref_cor[chimera,],spar = 0.8)
    smooth_spline_ko = smooth.spline(x = dev_time,y = ko_ref_cor[chimera,],spar = 0.8)
    
    png(sprintf("%s/time_dist_cor_%s.png",fig_dir,chimera))
    plot(x = dev_time,y= host_ref_cor[chimera,],pch = 19,main = chimera,ylim = c(min(host_ref_cor,ko_ref_cor),1),
         xlab = "developmental time",cex = 0.7,ylab= "Correlation time distributions")
    points(x =dev_time,y = ko_ref_cor[chimera,],pch = 19,col = "coral2",cex = 0.7)
    lines(x = dev_time,y = smooth_spline_host$y)
    abline(v = dev_time[which.max(smooth_spline_host$y)],lty = "dashed")
    lines(x = dev_time,y = smooth_spline_ko$y,col = "coral2")
    abline(v = dev_time[which.max(smooth_spline_ko$y)],lty = "dashed",col = "coral2")
    legend(x = "topleft",legend = c("KO","host"),pch = 19,col = c("coral2","black"))
    dev.off()
    
    
  }

}



