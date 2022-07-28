library("dplyr")


merge_umi_mat_with_embexe = function(scmat) {
  
  mat1 = scmat
  mat2 = scdb_mat("embexe")
  
  md2 = mat2@cell_metadata
  md2$cell_type = ifelse(md2$embryo == "empty","empty","WT")

  ignored_cls = union(mat1@ignore_cells,mat2@ignore_cells)
  ignored_genes = union(mat1@ignore_genes,mat2@ignore_genes)
  
  mat_all = rbind(cbind(mat1@mat,mat1@ignore_cmat),cbind(mat1@ignore_gmat,mat1@ignore_gcmat))
  
  mat_tmp = rbind(cbind(mat2@mat,mat2@ignore_cmat),cbind(mat2@ignore_gmat,mat2@ignore_gcmat))
  mat_tmp = mat_tmp[rownames(mat_all),]
  mat_all = cbind(mat_all,mat_tmp)
  md_all = bind_rows(mat1@cell_metadata,md2)
  
  rownames(md_all) = md_all$cell
  
  mat_new = scm_new_matrix(mat = mat_all, stat_type =  "umi",cell_metadata = md_all)
  mat_new = scm_ignore_cells(scmat = mat_new,ig_cells = ignored_cls)
  mat_new = scm_ignore_genes(scmat = mat_new,ig_genes = ignored_genes)
  
  
  return(mat_new)
}