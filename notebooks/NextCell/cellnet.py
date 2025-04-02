import torch
import torch.nn as nn
import torch.optim as optim
from transformers import BertModel, BertTokenizer
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
import time
import scipy.sparse as sparse
from scipy.spatial.distance import pdist
import glob
import networkx as nx
from importlib import reload
from sklearn_ann.kneighbors.annoy import AnnoyTransformer
from sklearn.metrics import mean_squared_error

import anndata as an
import scanpy as sc
import scanpy.external as sce
import phate


class CellNet:
    """
    A class to encapsulate AnnData to NetworkX graph conversion and related operations.
    """

    def __init__(
        self, 
        adata, 
        label_column='cell_type', 
        graph_key='X_pca',
        graph_neighbors=5,
        graph_n_components=2, 
        graph_metric='euclidean',
        distance_key='X_pca',
        distance_metric='euclidean',
        distance_n_components=2,
        source_label=None,
        target_label=None,
        verbose=True
    ):
        """
        Initializes the CellNet object
        """
        self.adata = adata
        self.index = adata.obs_names
        self.label_column = label_column
        self.source_label = source_label
        self.target_label = target_label
        self.verbose = verbose
        
        self.graph_key = graph_key
        self.graph_neighbors = graph_neighbors
        self.graph_n_components = graph_n_components
        self.graph_metric = graph_metric
        self.graph = self.build_graph()

        self.distance_key = distance_key
        self.distance_metric = distance_metric
        self.distance_n_components = distance_n_components
        
    def build_graph(self):
        """
        Builds a graph from the AnnData object.
    
        Args:
            verbose (bool, optional): If True, prints timing information. Defaults to False.
        """
        start_time = time.time()
    
        sc.pp.neighbors(
            self.adata,
            n_neighbors=self.graph_neighbors,
            use_rep=self.graph_key,
            n_pcs=self.graph_n_components,
            metric=self.graph_metric,
        )
    
        neighbors_time = time.time()
        if self.verbose:
            print(f"Neighbors calculation time: {neighbors_time - start_time:.4f} seconds")
    
        G = nx.from_scipy_sparse_array(
            self.adata.obsp['connectivities'],
        )
    
        graph_time = time.time()
        if self.verbose:
            print(f"Graph creation time: {graph_time - neighbors_time:.4f} seconds")
    
        # relabel nodes
        node_mapping = {i: self.index[i] for i in range(len(self.index))}
        G = nx.relabel_nodes(G, node_mapping)
        
        if self.verbose:
            print(G)
            print(f"Graph built.\n")

        return G

    def _random_source(self):
        if self.source_label is None:
            raise ValueError('CellNet.source_label cannot be `None`. Please set a source_label.')
            
        return random.choice(
            self.adata.obs[self.adata.obs[self.label_column] == self.source_label].index
        )

    def _random_target(self):
        if self.target_label is None:
            raise ValueError('CellNet.target_label cannot be `None`. Please set a target_label.')
            
        return random.choice(
            self.adata.obs[self.adata.obs[self.label_column] == self.target_label].index
        )
    

    def get_random_nodes(self, source_type, target_type):
        """
        Selects a random source node and a random target node from the AnnData object.
        """
        source_node = random.choice(self.adata.obs[self.adata.obs[self.annotation_col] == source_type].index)
        target_node = random.choice(self.adata.obs[self.adata.obs[self.annotation_col] == target_type].index)
        return source_node, target_node
    
    def _distance(self, node_i, node_j):
        """
        Computes the hueristic distance between two nodes.
        """
        i = self.index.get_loc(node_i)
        j = self.index.get_loc(node_j)
        key = self.distance_key
        X = self.adata[[i, j], :].obsm[key]
        return float(pdist(X, metric=self.distance_metric)[0])


    def find_path(self, source_node=None, target_node=None):
        path = nx.astar_path(
            self.graph,
            source_node,
            target_node,
            heuristic=lambda u, v: self._distance(u, v)
        )
        return path
    
    def get_random_paths(self, n_paths=1):
        """A function to generate n random paths between source and target types,
        with optional timing print statements.
        """
    
        paths = []
    
        for i in range(n_paths):
            start_time = time.time()  
            path = self.find_path(
                self._random_source(),
                self._random_target(),
            )
            end_time = time.time()  
            elapsed_time = end_time - start_time
    
            if self.verbose:
                print(f"Path {i+1} ({len(path)} nodes) generated in {elapsed_time:.4f} seconds.")
    
            paths.append(path) 
        return paths
  