import os 
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
# from torchvision import datasets, transforms
from StateNet import StateNet
from data_loader import Dataset

from data_loader import load_array
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from Config import Config, DataConfig

from utils import *

Config = Config()
DataConfig = DataConfig()

USE_CUDA = True if torch.cuda.is_available() else False
BATCH_SIZE = DataConfig.batch_size
N_EPOCHS = Config.epochs


from d2l import torch as d2l

from torch.utils import data
from data_loader import MyData

import torch
from sklearn.model_selection import KFold


def train(model, optimizer, train_loader, test_loader, epoch, classes_num=10):
    metric = d2l.Accumulator(3)
    capsule_net = model
    capsule_net.train()
    n_batch = len(list(enumerate(train_loader)))
    for batch_id, (data, target) in enumerate(train_loader):
        target = torch.sparse.torch.eye(classes_num).index_select(dim=0, index=target)

        if USE_CUDA:
            data, target = data.cuda(), target.cuda()
        
        X_ls = X2X_ls(data)

        optimizer.zero_grad()
        output, reconstructions, masked, _ = capsule_net(X_ls)
        loss = capsule_net.loss(data, output, target, reconstructions)
        loss.backward()
        optimizer.step()
        correct = sum(np.argmax(masked.cpu().numpy(), 1) == np.argmax(target.cpu().numpy(), 1))
        train_loss = loss.item()
        with torch.no_grad():
            metric.add(train_loss * data.shape[0], correct, data.shape[0])
        train_l = metric[0] / metric[2]
        train_acc = metric[1] / metric[2]
        if (batch_id + 1) % (n_batch // 5) == 0 or batch_id == n_batch - 1:
            animator.add(epoch + (batch_id + 1) / n_batch,
                            (train_l, train_acc, None))
    
    test_l, test_acc = test(capsule_net, test_loader, epoch, classes_num=Config.classes)

    print("Training at Epoch " + str(epoch))
    if epoch > 180:
        with_reconstruction=True
        torch.save(capsule_net.state_dict(),
                'results/{:03d}_model_dict_{}routing_reconstruction{}.pth'.format(epoch, Config.routing_iterations, with_reconstruction))


    return train_acc, train_l, test_l, test_acc

def test(capsule_net, test_loader, epoch, classes_num=10):
    capsule_net.eval()
    metric = d2l.Accumulator(3)
    with torch.no_grad():
        for batch_id, (data, target) in enumerate(test_loader):

            target = torch.sparse.torch.eye(classes_num).index_select(dim=0, index=target)
            
            if USE_CUDA:
                data, target = data.cuda(), target.cuda()
            
            X_ls = X2X_ls(data)

            output, reconstructions, masked, c_ij = capsule_net(X_ls)
            loss = capsule_net.loss(data, output, target, reconstructions)
            test_loss = loss.item()
            correct = sum(np.argmax(masked.cpu().numpy(), 1) ==
                            np.argmax(target.cpu().numpy(), 1))
            
            metric.add(test_loss, correct, data.shape[0])
            test_l = metric[0] / metric[2]
            test_acc = metric[1] / metric[2]

        animator.add(epoch + 1, (None, None, test_acc))
    return test_l, test_acc




capsule_net = StateNet(config=Config)

capsule_net = torch.nn.DataParallel(capsule_net)

if USE_CUDA:
    capsule_net = capsule_net.cuda()
capsule_net = capsule_net.module

optimizer = torch.optim.Adam(
        capsule_net.parameters(), 
        lr = Config.lr
    )


x_train, x_test, y_train, y_test = load_array(DataConfig)

# 5fold
k_folds = 5
kf = KFold(n_splits=k_folds, shuffle=True)

for fold, (train_idx, val_idx) in enumerate(kf.split(x_train)):
    x_train_fold = x_train[train_idx]
    y_train_fold = y_train[train_idx]
    x_val_fold = x_train[val_idx]
    y_val_fold = y_train[val_idx]
    train_dataset = MyData(x_train_fold, y_train_fold)
    train_loader = data.DataLoader(train_dataset, DataConfig.batch_size, shuffle=False)
    val_dataset = MyData(x_val_fold, y_val_fold)
    val_loader = data.DataLoader(val_dataset, DataConfig.batch_size, shuffle=False)

    
    train_acc_ls = []
    train_loss_ls = []
    val_acc_ls = []
    val_loss_ls = []

    animator = d2l.Animator(xlabel='epoch',xlim=[1, Config.epochs],
                    legend=['train loss', 'train acc', 'val acc'])
    for epoch in range(1, N_EPOCHS + 1):
        train_acc, train_l, val_l, val_acc = train(capsule_net, optimizer, train_loader, val_loader, epoch, classes_num=Config.classes)
        print(train_acc)
        print(train_l)
        train_acc_ls.append(train_acc)
        train_loss_ls.append(train_l)
        val_acc_ls.append(val_acc)
        val_loss_ls.append(val_l)

    print(f'loss {train_l:.3f}, train acc {train_acc:.3f}, '
        f'val acc {val_acc:.3f}')

    with open("./results/"+"Fold"+str(fold+1)+"_train_acc_ls.txt", "w") as f:
        for i in train_acc_ls:
            f.write(str(i)+'\n')
    with open("./results/"+"Fold"+str(fold+1)+"_train_loss_ls.txt", "w") as f:
        for i in train_loss_ls:
            f.write(str(i)+'\n')
    with open("./results/"+"Fold"+str(fold+1)+"_val_acc_ls.txt", "w") as f:
        for i in val_acc_ls:
            f.write(str(i)+'\n')
    with open("./results/"+"Fold"+str(fold+1)+"_val_loss_ls.txt", "w") as f:
        for i in val_loss_ls:
            f.write(str(i)+'\n')
            
    print(f'Fold {fold+1}: Validation Loss = {val_l:.3f}, Validation Accuracy = {val_acc:.3f}')
    torch.save(capsule_net.state_dict(), f'./results/{fold+1}_model_dict.pth')
    torch.save(capsule_net(), f'./results/{fold+1}_model.pth')
    plt.savefig('./results/'+'Fold'+str(fold+1)+'_Train_plot.png', dpi=300)
    plt.savefig('./results/'+'Fold'+str(fold+1)+'_Train_plot.pdf', dpi=300)
    optimizer.zero_grad()



test_dataset = MyData(x_test, y_test) 
test_loader = data.DataLoader(test_dataset, DataConfig.batch_size, shuffle=False)

capsule_net.eval()
test_l, test_acc = test(capsule_net, test_loader, 0, classes_num=Config.classes)

with open("./results/Test_results.txt", 'w') as f:
    f.write(f'TEST LOSS {test_l:.3f}, TEST ACCURACY {test_acc:.3f}')

