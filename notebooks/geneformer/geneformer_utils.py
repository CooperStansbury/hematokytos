import pandas as pd 
import numpy as np
import anndata as an
import scanpy as sc
import pickle
import torch

from datasets import Dataset, load_from_disk, load_dataset
import geneformer

def extract_gene_embeddings_batched(model, data, token_gene_dict, layer_to_quant=-1, process_batch=500, gpu_batch=10, verbose=False):
    """Extracts gene embeddings from a dataset in batches using a provided model.
    
    Args:
        model: The pre-trained model used for embedding generation.
        data: The input dataset (datasets.arrow_dataset.Dataset).
        token_gene_dict: A dictionary mapping token IDs to gene IDs.
        layer_to_quant: The model layer to extract embeddings from (default: -1, the last layer).
        process_batch: The number of samples to process in each outer loop iteration.
        gpu_batch: The batch size for the inner embedding generation (on the GPU).
    
    Returns:
        An AnnData object containing the gene embeddings and associated metadata.
    """
    total_batch_length = len(data)

    if verbose:
        print(f"Input data shape {data.shape}")
    
    all_embs = []
    for i in range(0, total_batch_length, process_batch):
        if verbose:
            print(f"Working batch {i+1}...")
        max_range = min(i + process_batch, total_batch_length)
        minibatch = data.select([i for i in range(i, max_range)])
        
        embs = geneformer.emb_extractor.get_embs(
            model=model,
            filtered_input_data=minibatch,  
            emb_mode='gene',
            layer_to_quant=-1,
            pad_token_id=0, 
            forward_batch_size=gpu_batch, 
            token_gene_dict=token_gene_dict, 
            special_token=False,
            summary_stat=None,
            silent=True,  
        )
        embs = embs.cpu().numpy()
        torch.cuda.empty_cache()
        all_embs.append(embs)

    
    all_embs = np.concatenate(all_embs, axis=0) # concatenate batches
    if verbose:
        print(f"{all_embs.shape=}")
    
    # # structure obs and prepare X
    X = []
    obs = data.to_pandas().reset_index(names='cell_id')
    if verbose:
        print(f"Raw obs: {obs.shape}")

    for i, row in obs.iterrows():
        length = row['length']
        cell_embedding = all_embs[i, :length, :]
        X.append(cell_embedding)
    
    X = np.concatenate(X, axis=0)
    obs = obs.explode('input_ids').reset_index(drop=True)
    obs['gene_id'] = obs['input_ids'].map(token_gene_dict)
    if verbose:
        print(f"Processed obs: {obs.shape}")
        print(f"Processed X: {X.shape}")
    edata = an.AnnData(X, obs=obs)
    return edata


def extract_cell_embeddings(model, data, token_gene_dict, layer_to_quant=-1, forward_batch_size=10):
    """
    Extracts cell embeddings from a model and returns them as an AnnData object.
    
    Args:
        model: The model to use for embedding extraction.
        data: The input data.
        token_gene_dict: Maps tokens to gene identifiers.
        layer_to_quant: The layer to extract embeddings from (default: last layer).
        forward_batch_size: Batch size for forward passes (default: 10).
    
    Returns:
        anndata.AnnData: An AnnData object containing the extracted cell embeddings.
    """
    embs = geneformer.emb_extractor.get_embs(
        model=model,
        filtered_input_data=data,
        emb_mode='cell',
        layer_to_quant=layer_to_quant,
        pad_token_id=0,  # Assuming this is a constant parameter for the function
        forward_batch_size=forward_batch_size,
        token_gene_dict=token_gene_dict,
        special_token=False,
        summary_stat=None,  
        silent=False, 
    )
    edata = an.AnnData(embs.cpu().numpy())
    edata.obs = data.to_pandas().astype(str)
    return edata


def extract_gene_embeddings(model, data, token_gene_dict, layer_to_quant=-1, forward_batch_size=10):
    """Extracts embeddings from a model and returns them as a DataFrame.

    This function provides an in-memory extraction of embeddings, allowing for convenient
    manipulation and analysis directly within your Python environment.

    Args:
        model: The model to use for embedding extraction.
        data: The input data for which embeddings need to be extracted.
        token_gene_dict: the token dictionary
        layer_to_quant (int, optional): The layer to quantize. Defaults to -1 (last layer).
        forward_batch_size (int, optional): The batch size for forward passes. Defaults to 10.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted embeddings.

    Raises:
        TypeError: If `model` is not a supported model type.
        ValueError: If `data` is not in the correct format.
    """
    embs = geneformer.emb_extractor.get_embs(
        model=model,
        filtered_input_data=data,
        emb_mode='gene',
        layer_to_quant=layer_to_quant,
        pad_token_id=0,  # Assuming this is a constant parameter for the function
        forward_batch_size=forward_batch_size,
        token_gene_dict=token_gene_dict,
        special_token=False,
        summary_stat=None,  
        silent=False, 
    )
    embs = embs.cpu().numpy()

    # prepare X
    n_cells, n_genes, n_dims = embs.shape
    X = embs.reshape(n_cells * n_genes, n_dims)

    # structure obs
    obs = data.to_pandas()
    obs = obs.explode('input_ids').reset_index(drop=True)
    obs['gene_id'] = obs['input_ids'].map(token_gene_dict)
    edata = an.AnnData(X, obs=obs.astype(str))
    return edata



def load_pickle(path):
    """Loads a pickled object from the specified file path.

    Args:
        path (str): The file path to the pickled object.

    Returns:
        The unpickled object.

    Raises:
        FileNotFoundError: If the file does not exist.
        pickle.UnpicklingError: If there's an error unpickling the object.
    """
    
    with open(path, "rb") as f:
        return pickle.load(f)

    
def load_model(model_path, model_type='Pretrained', n_classes=0, mode='eval'):
    """
    Loads a pre-trained or custom model for geneformer perturbations.

    Args:
        model_path (str): Path to the model file.
        model_type (str, optional): Type of model ('Pretrained' or custom). Default: 'Pretrained'.
        n_classes (int, optional): Number of output classes for custom models. Default: 0.
        mode (str, optional): Mode to load the model in ('eval' or 'train'). Default: 'eval'.

    Returns:
        The loaded model object.
    """

    model = geneformer.perturber_utils.load_model(
        model_type,
        n_classes,
        model_path,
        mode
    )

    return model


def load_data_as_dataframe(data_path, num_cells=None, shuffle=False) -> pd.DataFrame:
    """Loads a dataset, optionally shuffles it, and returns a subset as a Pandas DataFrame.

    Args:
        data_path: Path to the dataset file.
        num_cells: Number of cells to include in the subset (default: 100).
        shuffle: Whether to shuffle the dataset before subsetting (default: True).

    Raises:
        ValueError: If the requested subset size exceeds the dataset length.

    Returns:
        The subset of data as a Pandas DataFrame.
    """

    data = load_from_disk(data_path)

    if shuffle:
        data = data.shuffle(seed=42)

    if num_cells is None:
        return data.to_pandas()
    elif num_cells > len(data):
        raise ValueError(f"Requested subset size ({num_cells}) exceeds dataset length ({len(data)}). For all cells, use num_cells=`None.'")
    else:
        data_subset = data.select([i for i in range(num_cells)])
        return data_subset.to_pandas()
    

def make_embedding_anndata(embedding_df, data):
    """A function to make an anndata object of embeddings"""

    adata = an.AnnData(embedding_df.to_numpy())
    adata.obs = data
    return adata
    
    
def embedding_to_adata(df: pd.DataFrame, n_dim: int = None) -> an.AnnData:
    """Converts a Pandas DataFrame with an embedding to an AnnData object.

    Args:
        df: The input DataFrame with numerical embedding columns and optional metadata columns.
        n_dim: The number of dimensions to keep in the embedding. If None, all dimensions are kept.

    Returns:
        The converted AnnData object.

    Raises:
        ValueError: If `n_dim` exceeds the available dimensions in the DataFrame.
    """

    if n_dim is not None and n_dim > df.shape[1]:
        raise ValueError(f"n_dim ({n_dim}) exceeds available dimensions ({df.shape[1]})")

    # Assuming embedding columns are those that are not integers
    is_metadata = df.columns.astype(str).str.isdigit()
    metadata_df = df.loc[:, ~is_metadata]
    embedding_df = df.loc[:, is_metadata]

    cell_index = pd.Index([f"C{x}" for x in range(df.shape[0])], name='obs_names')

    if n_dim is not None:
        embedding_df = embedding_df.iloc[:, :n_dim]

    var_index = pd.Index([f"D{x}" for x in range(embedding_df.shape[1])], name='var_names')

    adata = an.AnnData(embedding_df.to_numpy())
    adata.obs_names = cell_index
    adata.var_names = var_index
    adata.obs = metadata_df
    return adata