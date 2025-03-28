# import torch
# import torch.nn as nn
# import torch.optim as optim
# from transformers import BertModel, BertTokenizer
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import random
# import glob
# import networkx as nx

# import anndata as an
# import scanpy as sc


# class NextCell(nn.Module):
#     def __init__(self, input_size=1000, hidden_size=128, output_size=1000, num_layers=2, num_heads=2, learning_rate=0.001):
#         super(NextCell, self).__init__()
#         self.input_size = input_size
#         self.hidden_size = hidden_size
#         self.output_size = output_size

#         self.transformer_encoder = nn.TransformerEncoder(
#             nn.TransformerEncoderLayer(input_size, num_heads, hidden_size),
#             num_layers
#         )
#         self.fc1 = nn.Linear(input_size, hidden_size)
#         self.relu = nn.ReLU()
#         self.fc2 = nn.Linear(hidden_size, output_size)

#         self.criterion = nn.MSELoss()
#         self.optimizer = optim.Adam(self.parameters(), lr=learning_rate)

#     def forward(self, input_gene_expression):
#         # Input shape: (batch_size, seq_len, input_size)
#         x = self.transformer_encoder(input_gene_expression.permute(1, 0, 2))  # (seq_len, batch_size, input_size)
#         x = x.permute(1, 0, 2).mean(dim=1)  # (batch_size, input_size)
#         x = self.fc1(x)
#         x = self.relu(x)
#         x = self.fc2(x)
#         return x

#     def train_model(self, input_gene_expression, target_gene_expression, epochs=10):
#         training_loss = []
#         for epoch in range(epochs):
#             self.optimizer.zero_grad()
#             predicted_gene_expression = self.forward(input_gene_expression)
#             loss = self.criterion(predicted_gene_expression, target_gene_expression)
#             loss.backward()
#             self.optimizer.step()
#             training_loss.append({'epoch': epoch + 1, 'loss': loss.item()})
#             print(f'Epoch {epoch + 1}, Loss: {loss.item()}')
#         return pd.DataFrame(training_loss)

#     def predict(self, input_gene_expression):
#         self.eval()
#         with torch.no_grad():
#             predictions = self.forward(input_gene_expression)
#         self.train()
#         return predictions

#     def get_embeddings(self, input_gene_expression):
#         self.eval()
#         with torch.no_grad():
#             embeddings = self.transformer_encoder(input_gene_expression.permute(1, 0, 2))
#             embeddings = embeddings.permute(1, 0, 2).mean(dim=1)
#         self.train()
#         return embeddings