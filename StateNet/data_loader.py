import pandas as pd
import numpy as np
import random
from torch.utils import data
from torch.utils.data import Dataset
from sklearn.model_selection import train_test_split

from Config import Config, DataConfig

DataConfig = DataConfig()
Config = Config()





class MyData(Dataset):
    def __init__(self, features, labels):
        self.features = features
        self.labels = labels

    def __getitem__(self, index):
        return self.features[index], self.labels[index]
    
    def __len__(self):
        return len(self.features)


def load_array(DataConfig):

    features = pd.read_csv(DataConfig.data_path, dtype=np.float32)
    labels = pd.read_csv(DataConfig.label_path)
    labels_oh = pd.read_csv(DataConfig.onehot_path)

    testing = True
    # testing = False
    if  testing:
        subset_size = 1000
        sp_idx = random.sample(range(labels.shape[0]), subset_size)
        features = features.iloc[sp_idx,]
        labels = labels.iloc[sp_idx,]
    
    features = np.array(features, dtype=np.float32)
    labels = pd.merge(labels, labels_oh, how='left', on=None, left_on=None, right_on=None).iloc[:,1]
    labels = np.array(labels, dtype=np.int32)
    
    x_train, x_test, y_train, y_test = train_test_split(features, labels, test_size = 0.1, random_state= 123)
    return x_train, x_test, y_train, y_test


if __name__ == '__main__':

    DataConfig = DataConfig()

    x_train, x_test, y_train, y_test = load_array(DataConfig)

    train_dataset = MyData(x_train, y_train)
    test_dataset = MyData(x_test, y_test)
    train_loader = data.DataLoader(train_dataset, DataConfig.batch_size, shuffle=False)
    test_loader = data.DataLoader(test_dataset, DataConfig.batch_size, shuffle=False)
    
    print(next(iter(train_loader)))




