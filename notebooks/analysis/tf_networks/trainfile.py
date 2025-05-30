import sys

sys.path.append('../')

from gears.utils import dataverse_download
from zipfile import ZipFile 
from gears import PertData, GEARS

# download pretrained model
## Download dataloader from dataverse
dataverse_download('https://dataverse.harvard.edu/api/access/datafile/6979957', 'norman_umi_go.tar.gz')

## Extract and set up dataloader directory
import tarfile
with tarfile.open('norman_umi_go.tar.gz', 'r:gz') as tar:
    tar.extractall()

## Download model from dataverse
dataverse_download('https://dataverse.harvard.edu/api/access/datafile/10457098', 'model.zip')

## Extract and set up model directory
with ZipFile(('model.zip'), 'r') as zip:
    zip.extractall(path = './')



# data_path = './hyb/'
# data_name = 'hybrid'
model_name = 'gears_misc_umi_no_test'

# pert_data = PertData(data_path)
# pert_data.load(data_path = data_path + data_name)
# pert_data.prepare_split(split = 'no_test', seed = 1)
# pert_data.get_dataloader(batch_size = 32, test_batch_size = 128)

import sys
sys.path.append('../')

from gears import PertData

# pert_data = PertData('./hyb_merged') # specific saved folder
# # pert_data.new_data_process(dataset_name = 'hybrid_norman', adata = adata_merged) # specific dataset name and adata object
# pert_data.load(data_path = './hyb_merged/hybrid_norman') # load the processed data, the path is saved folder + dataset_name
# pert_data.prepare_split(split = 'simulation', seed = 1) # get data split with seed
# pert_data.get_dataloader(batch_size = 32, test_batch_size = 128) # prepare data loader



# pert_data = PertData('./hyb_two') # specific saved folder
# # pert_data.new_data_process(dataset_name = 'hybrid_norman', adata = adata_merged) # specific dataset name and adata object
# pert_data.load(data_path = './hyb_two/hybridtwo') # load the processed data, the path is saved folder + dataset_name
# pert_data.prepare_split(split = 'simulation', seed = 1) # get data split with seed
# pert_data.get_dataloader(batch_size = 32, test_batch_size = 128) # prepare data loader


pert_data = PertData('./hyb_merged') # specific saved folder
# pert_data.new_data_process(dataset_name = 'hybrid_norman', adata = adata_merged) # specific dataset name and adata object
pert_data.load(data_path = './hyb_merged/hybrid_norman') # load the processed data, the path is saved folder + dataset_name
pert_data.prepare_split(split = 'simulation', seed = 1) # get data split with seed
pert_data.get_dataloader(batch_size = 32, test_batch_size = 128) # prepare data loader




gears_model = GEARS(pert_data, device = 'cuda', 
                        weight_bias_track = False, 
                        proj_name = 'pertnet', 
                        exp_name = 'pertnet')
gears_model.model_initialize(hidden_size = 64)

# gears_model.tunable_parameters()

# gears_model.train(epochs = 1, lr = 1e-3)

# gears_model.save_model('test_model')
gears_model.load_pretrained('test_model')


# gears_prediction = gears_model.predict([['FEV'], ['FEV', 'AHR']])
gears_prediction = gears_model.predict([['PRRX1'], ['PRRX1', 'MYOD1'], ['MYOD1']])


print(gears_prediction)
gears_prediction['genes'] = gears_model.gene_list


import pickle

with open('outputs/perturbation_full_with_genes.pickle', 'wb') as handle:
    pickle.dump(gears_prediction, handle, protocol=pickle.HIGHEST_PROTOCOL)

# gears_model = GEARS(pert_data, device = 'cuda', 
#                         weight_bias_track = False, 
#                         proj_name = 'gears', 
#                         exp_name = model_name)
# gears_model.load_pretrained('./model_ckpt')
# prediction = gears_model.predict([['CNN1', 'CBL']]) 

# print(f"printing prediction : {prediction}")


# ## If reproducing results from paper, you can use the same gene set, 
# ## although the function works even if GI_genes_file is set to None

# dataverse_download('https://dataverse.harvard.edu/api/access/datafile/6979958', 
#                    'genes_with_hi_mean.npy')

# prediction = gears_model.GI_predict(['CNN1', 'CBL'], GI_genes_file='./genes_with_hi_mean.npy')
# print(f"printing prediction : {prediction}")

# pert_data = PertData('./data')
# pert_data.load(data_name = 'norman')
# pert_data.prepare_split(split = 'simulation', seed = 1)
# pert_data.get_dataloader(batch_size = 32, test_batch_size = 128)


# gears_model = GEARS(pert_data, device = 'cuda', 
#                         weight_bias_track = False, 
#                         proj_name = 'pertnet', 
#                         exp_name = 'pertnet')
# gears_model.model_initialize(hidden_size = 64)


# gears_model.tunable_parameters()


# gears_model.train(epochs = 1, lr = 1e-3)


# gears_model.save_model('test_model')
# gears_model.load_pretrained('test_model')





# gears_model.predict([['FEV'], ['FEV', 'AHR']])

# gears_model.gene_list[:5]