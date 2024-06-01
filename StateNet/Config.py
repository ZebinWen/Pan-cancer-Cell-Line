
class Config:
    def __init__(self, routing_iterations=3, **kwargs):
        # Train
        self.epochs = 20
        self.lr = 0.01

        # Cell State Capsules
        self.dim_NMF_Program = 50
        self.pc_dim_caps = 32
        self.pc_num_caps = 21

        # Cell Line Capsules
        self.classes = 245
        self.dc_dim_caps = 16

        # routing by agreement
        self.routing_iterations = 1

        # # fluctuation operation
        # self.fluc=0.1
        # self.propick=0


class DataConfig:
    def __init__(self):
        self.batch_size = 128
        self.data_path = './data/1.InputForStateNet_CL.MP.total.mtx_ColGene_RowSample.csv'
        self.label_path = './data/1.InputLabel_label.ls.csv'
        self.onehot_path = './data/1.LabelOneHot.df.csv'



