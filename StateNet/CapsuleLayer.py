from utils import *
import torch
from torch import nn
from torch.nn import functional as F
from torch.autograd import Variable
from torch.nn.parameter import Parameter

USE_CUDA = True if torch.cuda.is_available() else False


class StateCaps(nn.Module):
    def __init__(self, dim_NMF_Program=50, out_channels=8, pc_num_caps=21):
        super(StateCaps, self).__init__()

        self.capsules = nn.ModuleList([
            nn.Linear(in_features=dim_NMF_Program, out_features=out_channels)
            for _ in range(pc_num_caps)])

    def forward(self, x):
        u = [capsule(x[i]) for i, capsule in enumerate(self.capsules)]
        u = torch.stack(u, dim=1)
        return F.sigmoid(u)


class CellLineCaps(nn.Module):
    def __init__(self, fluc, propick, num_capsules=10, num_routes=32 * 6 * 6, in_channels=8, out_channels=16, num_iterations=3):
        super(CellLineCaps, self).__init__()
        self.fluc = fluc
        self.propick = propick
        self.in_channels = in_channels
        self.num_routes = num_routes
        self.num_capsules = num_capsules
        self.num_iterations = num_iterations
        self.W = Parameter(nn.init.kaiming_uniform_(torch.empty(1, num_routes, num_capsules, out_channels, in_channels)))

    def forward(self, x):
        batch_size = x.size(0)
        x = torch.stack([x] * self.num_capsules, dim=2).unsqueeze(4)

        W = torch.cat([self.W] * batch_size, dim=0) 
        u_hat = torch.matmul(W, x) 

        b_ij = Variable(torch.zeros(1, self.num_routes, self.num_capsules, 1)) 
        if USE_CUDA:
            b_ij = b_ij.cuda()

        for iteration in range(self.num_iterations):
            c_ij = F.softmax(b_ij, dim=1)
            c_ij = torch.cat([c_ij] * batch_size, dim=0).unsqueeze(4)
            if self.fluc is not None and self.propick is not None and iteration == self.num_iterations - 1:
                c_ij[:,self.propick] = c_ij[:,self.propick] + self.fluc

            s_j = (c_ij * u_hat).sum(dim=1, keepdim=True)
            v_j = squash(s_j)

            if iteration < self.num_iterations - 1:
                a_ij = torch.matmul(u_hat.transpose(3, 4), torch.cat([v_j] * self.num_routes, dim=1))
                b_ij = b_ij + a_ij.squeeze(4).mean(dim=0, keepdim=True)

        return v_j.squeeze(1), c_ij
    
    def disturb(self, fluc, propick):
        self.fluc = fluc
        self.propick = propick
        print("Disturbtion of C value is:", fluc)
        print("The target program index is:", propick)


class Decoder(nn.Module):
    def __init__(self, output_genes=784, classes_num=None, dc_dim_caps=16):
        super(Decoder, self).__init__()
        self.classes_num = classes_num
        self.output_genes = output_genes
        self.reconstraction_layers = nn.Sequential(
            nn.Linear(dc_dim_caps * classes_num, 512),
            nn.ReLU(inplace=True),
            nn.Linear(512, 1024),
            nn.ReLU(inplace=True),
            nn.Linear(1024, self.output_genes),
            nn.Sigmoid()
        )

    def forward(self, x, data=None):
        classes = torch.sqrt((x ** 2).sum(2))
        classes = F.softmax(classes, dim=1)

        _, max_length_indices = classes.max(dim=1)
        masked = Variable(torch.sparse.torch.eye(self.classes_num))
        if USE_CUDA:
            masked = masked.cuda()
        masked = masked.index_select(dim=0, index=Variable(max_length_indices.squeeze(1).data))
        t = (x * masked[:, :, None, None]).view(x.size(0), -1)
        reconstructions = self.reconstraction_layers(t)
        reconstructions = reconstructions.view(-1, self.output_genes)
        return reconstructions, masked



