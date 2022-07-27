# Synchronization of Placenta and Embryonic Development Underlies Mouse Gastrulation


This repository is the accompanying code for our paper on the development of extraembryonic ectoderm and embryonic tissues during mouse gastrulation. There is an [MCView shiny app](https://tanaylab.weizmann.ac.il/EmbExe) where you can interrogate the data. The code is splitted into jupyter notebooks that can be found in the notebook folder.

## Running the notebooks

Prior to any analysis, after cloning the repository, please download first the necessary data by running (in the root directory of the cloned repository):


```bash
R -e "source('scripts/download_data.R'); download_full_data()"
```

This will download all the necessary data including processed metacell objects necessary for generating the figures. The subfolder `data/umi.tables` contains all the UMI matrices of the scRNA-seq data used in the paper.

## Necessary R packages

The initialization script (`scripts/init.R`) loads automatically the necessary R packages to run the notebooks. The analysis was done using R 4.0.5 and the following packages:

- gridExtra 2.3
- qvalue 2.22.0
- ggrepel 0.9.1
- ggpubr 0.4.0   
- ggplot2 3.3.5
- zoo 1.8-9
- qlcMatrix 0.9.7
- sparsesvd 0.2  
- slam 0.1-49
- tidyr 1.2.0
- dplyr 1.0.9
- tgutil 0.1.13  
- tgstat 2.3.17
- metacell 0.3.7
- Matrix 1.3-4

## Notebook order 

Plots shown in the paper can be reproduced by running the notebooks.



### Analysis of wildtype ExE and embryonic manifold
For reproducing the analysis of the wildtype ExE and embryonic manifold you should run the following notebooks in that order.

1.  import_mars
2.  embexe_find_bad_genes
3.  embexe metacell construction
4.  embryo_temporal_ordering
5.  embexe_interpolate_time
6.  mc2d_projection_embexe
7.  split_embexe_into_emb_and_exe
8.  emb.estimation_of_proliferation_rates
9.  emb_generate_network
10. mc2d_projection_emb
11. mc2d_projection_exe

### Analysis of EXE-specific Elf5 KO embryos 

10. import_elf5
11. elf5_processing_summary

### Analysis of *ex utero* cultured embryos

12. import_10x_exutero
13. exutero_doublet_removal
14. exutero_f_find_bad_genes
15. exutero_f_generate_metacell


