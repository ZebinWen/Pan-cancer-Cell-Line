import torch
from torch.nn import functional as F
from torch import nn
from CapsuleLayer import StateCaps, CellLineCaps, Decoder

from Config import Config




class StateNet(nn.Module):
    def __init__(self, config=None, fluc=None, propick=None, *args, **kwargs):
        super(StateNet, self).__init__(*args, **kwargs)
        
        self.primary_capsules = StateCaps(
            dim_NMF_Program = config.dim_NMF_Program,
            out_channels = config.pc_dim_caps,
            pc_num_caps = config.pc_num_caps
        )
        self.digit_capsules = CellLineCaps(
            num_capsules = config.classes,
            num_routes = config.pc_num_caps,
            in_channels = config.pc_dim_caps,
            out_channels = config.dc_dim_caps,
            num_iterations = config.routing_iterations, 
            fluc=fluc, 
            propick=propick
        )
        self.decoder = Decoder(
            output_genes = config.pc_num_caps * config.dim_NMF_Program,
            classes_num = config.classes,
            dc_dim_caps = config.dc_dim_caps
        )
        self.mse_loss = nn.MSELoss()
    
    def forward(self, data):
        output, c_ij = self.digit_capsules(self.primary_capsules(data))
        reconstructions, masked = self.decoder(output, data)
        return output, reconstructions, masked, c_ij

    def loss(self, data, x, target, reconstructions):
        return self.margin_loss(x, target) + self.reconstruction_loss(data, reconstructions)

    def margin_loss(self, x, labels):
        batch_size = x.size(0)

        v_c = torch.sqrt((x ** 2).sum(dim=2, keepdim=True))

        left = F.relu(0.9 - v_c).view(batch_size, -1)
        right = F.relu(v_c - 0.1).view(batch_size, -1)

        loss = labels * left + 0.5 * (1.0 - labels) * right 
        loss = loss.sum(dim=1).mean()

        return loss
    
    def reconstruction_loss(self, data, reconstructions):
        loss = self.mse_loss(reconstructions.view(reconstructions.size(0), -1), data.view(reconstructions.size(0), -1))
        return loss * 0.0005
    



if __name__ == '__main__':

    USE_CUDA = True if torch.cuda.is_available() else False
    Config = Config()

    capsule_net = StateNet(config=Config)

    capsule_net = torch.nn.DataParallel(capsule_net)

    if USE_CUDA:
        capsule_net = capsule_net.cuda()
    capsule_net = capsule_net.module

    print(capsule_net)

