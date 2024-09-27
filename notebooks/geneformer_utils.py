import pandas as pd 
import numpy as np
import anndata as an
import scanpy as sc
import pickle

from datasets import Dataset, load_from_disk, load_dataset
import geneformer


DEFAULT_NAME_PATH = "/nfs/turbo/umms-indikar/shared/projects/geneformer/geneformer/gene_name_id_dict.pkl"
DEFAULT_TOKEN_PATH = "/nfs/turbo/umms-indikar/shared/projects/geneformer/token_dictionary.pkl"
DEFAULT_MEDIAN_PATH = "/nfs/turbo/umms-indikar/shared/projects/geneformer/geneformer/gene_median_dictionary.pkl"


def extract_embedding_in_mem(model, data, emb_mode='cell', layer_to_quant=-1, forward_batch_size=10):
    """Extracts embeddings from a model and returns them as a DataFrame.

    This function provides an in-memory extraction of embeddings, allowing for convenient
    manipulation and analysis directly within your Python environment.

    Args:
        model: The model to use for embedding extraction.
        data: The input data for which embeddings need to be extracted.
        emb_mode (str, optional): The embedding mode. Defaults to 'cell'.
        layer_to_quant (int, optional): The layer to quantize. Defaults to -1 (last layer).
        forward_batch_size (int, optional): The batch size for forward passes. Defaults to 10.

    Returns:
        pandas.DataFrame: A DataFrame containing the extracted embeddings.

    Raises:
        TypeError: If `model` is not a supported model type.
        ValueError: If `data` is not in the correct format.
    """

    embs = geneformer.emb_extractor.get_embs(
        model,
        data,
        emb_mode,
        layer_to_quant,
        0,  # Assuming this is a constant parameter for the function
        forward_batch_size,
        summary_stat=None,  
        silent=False, 
    )
    data = embs.cpu().numpy()
    if emb_mode=='cell':
        return pd.DataFrame(data)
    else:
        return data



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