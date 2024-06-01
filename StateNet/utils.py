import torch

def squash(tensor, dim=-1):
    squared_norm = (tensor ** 2).sum(dim=dim, keepdim=True)
    scale = squared_norm / (1 + squared_norm)
    epsilon = 1e-7
    return scale * tensor / torch.sqrt(squared_norm + epsilon)

def safe_norm(x, dim=-1, epsilon=1e-7):
    return ((x ** 2).sum(dim=dim) + epsilon) ** 0.5

def X2X_ls(X, sep=50):
    for i in range(0, X.shape[1], sep):
        globals()[f"X{int(i/sep) + 1}"] = torch.tensor(X[:,i:i+sep])
    X_ls = []
    for i in range(int(X.shape[1]/sep)):
        X_ls.append(globals()[f"X{i+1}"])
    
    return X_ls

