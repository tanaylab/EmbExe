library("metacell")
source("scripts/generic_mc.r")
tgconfig::override_params(config_file = "config/emb_exe.yaml",package = "metacell")

scdb_init("scrna_db",force_reinit = T)
scfigs_init("figs")

# first iteration without out filtered genes
# generate_mc(mat_nm, color_key=NA,recompute = T)
# then run embdyn_find_bad_genes.r and filter bad gene modules
# included bad genes in second iteration
bad_genes = read.table("data/embexe.bad_genes.txt",sep = "\t",stringsAsFactors = F)
bad_genes = bad_genes[,1]

mat_nm = "embexe"

generate_mc(mat_nm, color_key=NA,add_bad_genes = bad_genes,recompute = T)

mcatlas_annotate_mc_by_mc2mc_projection(atlas_id = "emb_exe",qmc_id = "embexe_bs500f",qmat_naming_type = "mars",new_qmc_id = "embexe_recolored",q_gset_id = "embexe",
                                        fig_cmp_dir = "figs/embexe_atlas_projection")
