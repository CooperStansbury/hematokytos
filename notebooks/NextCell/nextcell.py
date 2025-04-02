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

import anndata as an
import scanpy as sc

class NextCell(nn.Module):
    def __init__(self, input_size=1000, hidden_size=128, output_size=1000, num_layers=2, num_heads=2, learning_rate=0.001):
        super(NextCell, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size

        self.transformer_encoder = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(
                input_size, 
                num_heads, 
                hidden_size,
            ),
                num_layers
        )
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, output_size)

        self.criterion = nn.MSELoss()
        self.optimizer = optim.Adam(self.parameters(), lr=learning_rate)

    def forward(self, x):
        """
        Predicts the next cell expression vector.

        Args:
            x (torch.Tensor): Input cell expression vector of shape [batch_size, input_size].

        Returns:
            torch.Tensor: Predicted cell expression vector of shape [batch_size, output_size].
        """
        x = x.unsqueeze(1) # Shape [batch_size, 1, input_size]
        x = x.transpose(0, 1) # Shape [1, batch_size, input_size] for transformer
        encoded = self.transformer_encoder(x) # Shape [1, batch_size, input_size]
        encoded = encoded.transpose(0, 1) # Shape [batch_size, 1, input_size]
        encoded = encoded.squeeze(1) # Shape [batch_size, input_size]

        x = self.fc1(encoded)
        x = self.relu(x)
        x = self.fc2(x)
        return x

    def train_step(self, x, y):
        """
        Performs a single training step.

        Args:
            x (torch.Tensor): Input cell expression vector of shape [batch_size, input_size].
            y (torch.Tensor): Target cell expression vector of shape [batch_size, output_size].

        Returns:
            float: Loss value.
        """
        self.optimizer.zero_grad()
        outputs = self.forward(x)
        loss = self.criterion(outputs, y)
        loss.backward()
        self.optimizer.step()
        return loss.item()

    def train_model(self, input_gene_expression, target_gene_expression, epochs=10, verbose=False):
        """
        Trains the model for a given number of epochs.

        Args:
            input_gene_expression (torch.Tensor): Input gene expression data.
            target_gene_expression (torch.Tensor): Target gene expression data.
            epochs (int): Number of training epochs.

        Returns:
            pd.DataFrame: Training loss for each epoch.
        """
        training_loss = []
        for epoch in range(epochs):
            loss = self.train_step(input_gene_expression, target_gene_expression)
            training_loss.append({'epoch': epoch + 1, 'loss': loss})
            if verbose:
                print(f'Epoch {epoch + 1}, Loss: {loss}')
        return pd.DataFrame(training_loss)

    def predict(self, input_gene_expression):
        """
        Predicts the next cell expression vector.

        Args:
            input_gene_expression (torch.Tensor): Input gene expression data.

        Returns:
            torch.Tensor: Predicted gene expression data.
        """
        self.eval()
        with torch.no_grad():
            predictions = self.forward(input_gene_expression)
        self.train()
        return predictions

    def get_embeddings(self, input_gene_expression):
        """
        Gets the embeddings from the transformer encoder.

        Args:
            input_gene_expression (torch.Tensor): Input gene expression data.

        Returns:
            torch.Tensor: Embeddings.
        """
        self.eval()
        with torch.no_grad():
            x = input_gene_expression.unsqueeze(1)
            x = x.transpose(0, 1)
            embeddings = self.transformer_encoder(x)
            embeddings = embeddings.transpose(0, 1).mean(dim=1)
        self.train()
        return embeddings

