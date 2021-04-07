import os
import scanpy as sc
from ._plotting import pair_plot, plot_z_3d, keys_to_colors
import numpy as np
from ._mymodel import MyModel
from . import setup_anndata


def save_results_txt(
    model_path, latent, cell_type, start_idx, end_idx, time_course, init_indices
):
    np.savetxt(os.path.join(model_path, "latent.txt"), latent)
    np.savetxt(os.path.join(model_path, "cell_type.txt"), cell_type, fmt="%s")
    np.savetxt(os.path.join(model_path, "start_idx.txt"), start_idx, fmt="%d")
    np.savetxt(os.path.join(model_path, "end_idx.txt"), end_idx, fmt="%d")
    np.savetxt(os.path.join(model_path, "time_course.txt"), time_course)
    np.savetxt(os.path.join(model_path, "init_indices.txt"), init_indices, fmt="%d")

    # print useful command
    cmd = "cp -r " + model_path + " ../../diffeq_match/case_studies/"
    print("useful command:\n", cmd)


def init_model_name(idx, m_type, n_latent, n_hidden, min_counts, n_top_genes):
    model_dir = "models"  # hardcoded
    if not os.path.exists(model_dir):
        os.mkdir(model_dir)
    model_name = (
        m_type
        + str(idx)
        + "_D"
        + str(n_latent)
        + "_G"
        + str(n_top_genes)
        + "_H"
        + str(n_hidden)
        + "_MC"
        + str(min_counts)
    )
    return model_name


def parse_args(data_dir, idx, model_name, pretrained):

    # Command line arguments
    idx = int(idx)
    pretrained = True if (int(pretrained) > 0) else False
    print("idx =", idx)
    print("pretrained =", pretrained)
    print(" ")

    # Read data and get final model name
    all_names = sorted(os.listdir(data_dir))
    for j, name in enumerate(all_names):
        print(j + 1, ":", name)
    fn = all_names[idx - 1]
    data_name = fn.split(".")[0]
    model_name = model_name + "-" + data_name

    adata = sc.read_h5ad(os.path.join(data_dir, fn))
    print("\nMODEL_NAME: >>>>>>>>>>", model_name, "<<<<<<<<<<<")
    return model_name, adata, pretrained


def analysis(
    data_dir, idx, pretrained, linear, n_latent, n_hidden, min_counts, n_top_genes
):

    # Setup
    m_type = "lin" if linear else "nl"
    model_name = init_model_name(
        idx, m_type, n_latent, n_hidden, min_counts, n_top_genes
    )
    model_name, adata, pretrained = parse_args(data_dir, idx, model_name, pretrained)

    # Get names of start and end cell(s), and pseudotime
    cell_ids = adata.obs["cell_id"]
    cell_type = adata.obs["cell_type"]
    start_idx = np.where(np.isin(cell_ids, adata.uns["start_id"]))[0]
    end_idx = np.where(np.isin(cell_ids, adata.uns["end_id"]))[0]
    time_course = np.array(adata.obs["timecourse"])
    print("start cell(s):", start_idx)
    print("end cell(s):", end_idx)

    # Get indices of all start cells
    init_type = np.array(cell_type[start_idx])[0]
    print("init_type: ", init_type)
    init_indices = np.where(np.array(cell_type) == init_type)[0]
    print("found", len(init_indices), "cells of type", init_type)

    # Set output file names
    pairs_fn = "pairs-" + model_name + ".png"
    umap1_fn = "-celltype-" + model_name + ".png"
    umap2_fn = "-timecourse-" + model_name + ".png"

    # Filter genes
    sc.pp.filter_genes(adata, min_counts=min_counts)

    # Normalize (store also raw)
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # freeze the state in `.raw`

    # Select highly variable genes
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, subset=True, layer="counts", flavor="seurat_v3"
    )

    # Setup adata and model directory
    setup_anndata(adata, layer="counts", labels_key="cell_type")
    model_path = os.path.join("models", model_name)
    print(adata)
    if not pretrained:
        # Create and train model
        model = MyModel(adata, n_latent=n_latent, n_hidden=n_hidden)
        model.train()
        model.save(model_path)
    else:
        # Load pretrained model
        model = MyModel(adata, n_latent=n_latent, n_hidden=n_hidden)
        model.load(model_path, adata, use_gpu=False)
        print("Loaded model from:", model_path)

    # Get latent and plot
    latent = model.get_latent_representation()
    cell_colors = keys_to_colors(cell_type)
    pair_plot(
        latent, cell_colors, save_name=pairs_fn, start_idx=start_idx, end_idx=end_idx
    )
    if n_latent == 3:
        plot_z_3d(latent, cell_colors)

    # Rec error
    rerr = model.get_reconstruction_error()
    print("Reconstruction error = ", rerr)

    adata.obsm["X_scVI"] = latent

    # use scVI latent space for UMAP generation
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.2)

    sc.pl.umap(adata, color=["cell_type"], frameon=False, save=umap1_fn, show=False)

    sc.pl.umap(adata, color=["timecourse"], frameon=False, save=umap2_fn, show=False)

    # save latent representation as txt
    save_results_txt(
        model_path, latent, cell_type, start_idx, end_idx, time_course, init_indices
    )

    return adata, model
