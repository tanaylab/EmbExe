#library("umap")
tgconfig::override_params(config_file = "config/emb_exe.yaml",package = "metacell")
# generate mgraph and mc2d object

# generate mgraph using logistic distance

feat_gset = "emb"
mc_id = "emb"
#mgraph_id = paste0(mc_id,"_logist")
mgraph_id = "emb_logist"
logist_loc = 1
logist_scale = 0.2
logist_eps = 4e-5
max_d_fold = 2
tgconfig::set_param(param = "mcell_mgraph_max_confu_deg",value = 4,package = "metacell")


mc = scdb_mc(mc_id)
gset = scdb_gset(feat_gset)
feat_genes = names(gset@gene_set)

mgraph = mgraph_comp_logist(mc = mc, genes = feat_genes, loc = logist_loc,scale =  logist_scale, eps = logist_eps, max_d_fold = max_d_fold)

scdb_add_mgraph(id = mgraph_id,mgraph = tgMCManifGraph(mc_id = mc_id,mgraph = mgraph))
