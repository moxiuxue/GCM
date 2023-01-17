import desc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
#%matplotlib inline 
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

adata= desc.read_10X("./10X")
adata.obs_names_make_unique()
info=pd.read_csv("./10X/batch.tsv", header=None, sep='\t')
adata.obs['batch']=list(info[1])

adata.raw=adata.copy()
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
desc.log1p(adata)
#adata.raw=adata.copy()
adata.X
sc.pp.highly_variable_genes(adata, n_top_genes=2000,batch_key='batch',subset=True)
desc.scale_bygroup(adata,groupby='batch')
adata = desc.train(adata, dims=[adata.shape[1], 256, 64], tol=0.01,
        n_neighbors=10,batch_size=256,
        louvain_resolution=[0.4,0.6,0.8,1.2,1.4,1.6,1.8,2.0],save_dir=hh+"result",
        do_tsne=True, learning_rate=300,do_umap=True, num_Cores_tsne=10,use_GPU=False,save_encoder_weights=False)

file_n = "_gcm_plot.pdf"

#Umap plot
sc.pl.scatter(adata,basis="umap0.4", legend_loc='on data',color=['desc_0.4','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap0.6", legend_loc='on data',color=['desc_0.6','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap0.8", legend_loc='on data',color=['desc_0.8','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap1.2", legend_loc='on data',color=['desc_1.2','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap1.4", legend_loc='on data',color=['desc_1.4','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap1.6", legend_loc='on data',color=['desc_1.6','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap1.8", legend_loc='on data',color=['desc_1.8','batch'],save=file_n)
sc.pl.scatter(adata,basis="umap2.0", legend_loc='on data',color=['desc_2.0','batch'],save=file_n)

adata.write(save_dir+hh+"_desc.h5ad")
