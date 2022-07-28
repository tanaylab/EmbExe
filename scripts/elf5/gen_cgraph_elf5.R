library("metacell")
tgconfig::override_params(config_file = "config/emb_exe.yaml",package = "metacell")

#scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")

gen_cgraph = function(mat_nm,T_vm = 0.1,Knn = 100) {
  
  bad_genes = read.table("data/embexe.bad_genes.txt",sep = "\t",stringsAsFactors = F)
  bad_genes = bad_genes[,1]
  bad_genes = c(bad_genes,"Elf5")
  
  mcell_add_gene_stat(mat_nm, mat_nm, force=T)
  
  mcell_gset_filter_varmean(mat_nm, mat_nm, T_vm=T_vm, force_new=T)
  mcell_gset_filter_cov(mat_nm, mat_nm, T_tot=50, T_top3=3)
  
  gset = scdb_gset(mat_nm)
  nms = names(gset@gene_set)
  #bad gene that will be removed from list of genes that helps to mark metacell 
  bad_g = c(grep("^Rpl",nms,v=T),grep("^Gm",nms,v=T),grep("Rps",nms,v=T))
  
  bad_g = c(bad_g, bad_genes)
  
  gset_f = gset_new_restrict_nms(gset=gset, bad_g, inverse=T, "feat filt")
  scdb_add_gset(mat_nm, gset_f)
  
  mcell_add_cgraph_from_mat_bknn(mat_id=mat_nm, 
                                 gset_id = mat_nm, 
                                 graph_id=mat_nm,
                                 K=Knn,
                                 dsamp=T)
  
}
