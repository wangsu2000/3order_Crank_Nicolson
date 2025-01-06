#coding=utf-8
import numpy as np
import cupy as cp
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from PIL import Image
import cv2
from tqdm import trange
import SimpleITK as sitk
import numpy as np
import os

### 生成高斯核权重
def gaussian_kernel(size, sigma):
    kernel = np.fromfunction(
        lambda x, y: (1/(2*np.pi*sigma**2)) * np.exp(-((x - size//2)**2 + (y - size//2)**2) / (2*sigma**2)),
        (size, size)
    )
    return kernel / np.sum(kernel)

### 生成高斯核
def create_gauss_kernel(sigma,kernel_size):
    gaussian = gaussian_kernel(kernel_size, sigma)
    gaussian = gaussian / np.sum(gaussian)
    kernel_gaussian = torch.reshape(torch.tensor(gaussian, dtype=torch.float32),(1,1,kernel_size,kernel_size))
    Gaussian = nn.Conv2d(1, 1, kernel_size=3, stride=1, padding=kernel_size // 2,padding_mode='zeros', bias=False)
    Gaussian.weight = nn.Parameter(kernel_gaussian, requires_grad=False)
    Gaussian.cuda()
    return Gaussian

def read_single_dicom_by_sitk(filepath):
    image = sitk.ReadImage(filepath)
    image_data = sitk.GetArrayFromImage(image)
    return np.squeeze(np.float32(image_data),axis = 0)

def rescale0_1(x,MIN_B=-1024, MAX_B=3072):
    x[x<-1024]=-1024
    x[x>3072]=3072
    return(x - MIN_B) / (MAX_B - MIN_B)

### 海氏函数
def Heaviside(phi,epsilon = 1.0):
    return 1.0 / np.pi * torch.atan(phi / epsilon) + 0.5

### 狄拉克函数
def Dirac(phi,epsilon = 1.0):
    return 1.0 / np.pi * epsilon / (epsilon**2 + phi**2)

### 更新b的函数
def getLatentb(phi,image,K,c1,c2):
    H = Heaviside(phi)
    numerator = K(image*H*c1 + image*(1.0 - H)*c2)
    denomenator = K(H*c1*c1 + (1.0 - H)*c2*c2)
    return numerator / denomenator

### 更新c的函数
def getLatentc(phi,image,K,b):
    H = Heaviside(phi)
    c1 = torch.sum(K(b)*image*H) / torch.sum(K(b**2)*H)
    c2 = torch.sum(K(b)*image*(1.0 - H)) / torch.sum(K(b**2)*(1.0 - H))
    return c1,c2

### LIC模型运行过程中需要的超参数
params = {
        'epsilon': 1.0,
        'nu': 0.001*255.0*255.0,
        'mu': 5.0e-2,
        'eta': 5.0e-1,
        'regulareps': 1.0e-14,
        'maxepoch':300,
        'default_c1':5.0e-1,
        'default_c2':-5.0e-1,
        'niter':8
    }

### LIC模型运行过程中每一步需要计算的梯度
def getGrad(phi,image,K,params):
    ## 这里应当循环
    c1 = params['default_c1']
    c2 = params['default_c2']
    for k in range(params['niter']):
        b = getLatentb(phi,image,K,c1,c2)
        c1,c2 = getLatentc(phi,image,K,b)
        if k > 0 and torch.sum((b - last_b)**2) < 1.0e-8 and torch.abs(last_c[0] - c1) < 1.0e-8 and torch.abs(last_c[1] - c2) < 1.0e-8:
#                 print("内部迭代收敛!,step at iter = {:d}".format(k+1))
                break
        last_b = b
        last_c = [c1,c2]
    ## 这里应当循环到收敛
    dirac = Dirac(phi)
    grad_image = ( -2*K(b)*image*c1 + K(b**2)*c1*c1 +\
                2*K(b)*image*c2 -K(b**2)*c2*c2  )*dirac
    px = torch.zeros_like(phi)
    px[:,:,1:-1,:] = (phi[:,:,2:,:] - phi[:,:,:-2,:]) / 2.0
    py = torch.zeros_like(phi)
    py[:,:,:,1:-1] = (phi[:,:,:,2:] - phi[:,:,:,:-2]) / 2.0
    denomenator = torch.sqrt(px**2 + py**2 + params['regulareps']) #加上这一项防止分母为0
    px = px / denomenator
    pxx = torch.zeros_like(phi)
    pxx[:,:,1:-1,:] = (px[:,:,2:,:] - px[:,:,:-2,:]) / 2.0
    py = py / denomenator
    pyy = torch.zeros_like(phi)
    pyy[:,:,:,1:-1] = (py[:,:,:,2:] - py[:,:,:,:-2]) / 2.0
    
    ### 得到平均曲率
    mean_curvature = pxx + pyy
    
    ### 计算的梯度没有包含Laplace项,便于半隐格式求解
    
    grad = grad_image - dirac*mean_curvature*params['nu'] - params['mu']*mean_curvature
    return grad

### 水平集初始化
def level_set_initialize(image,radius = 5):
    # 输入参数中的radius为mask为白色区域的半径,随机取中点
    # image 应当为[batchsize,1,width,height] 形状的
    volume = []
#     print(image.shape)
    xx,yy = np.meshgrid(np.arange(0,image.shape[-2]),np.arange(0,image.shape[-1]))
    for k in range(image.shape[0]):
        center_x = np.random.randint(low=0, high=image.shape[-2])
        center_y = np.random.randint(low=0, high=image.shape[-1])
        zz = (xx - center_x)**2 + (yy - center_y)**2 <= radius**2
        zz = zz.astype(np.float32)*2 - 1
        volume.append(zz.T)
    volume = np.stack(volume,axis = 0)
#     print(volume.shape)
    volume = torch.tensor(volume,dtype = torch.float32).unsqueeze(1).cuda()
#     print(volume.shape)
    return volume

### 用于梯度下降的轮子:
def fft_batch_helmholtz(b,lamb,h = 1.0):
        bs, m, n = b.shape
        row = n * cp.cos(cp.pi * cp.arange(0, n) / n)
        row[0] = row[0] + n
        col = m / 2 * cp.ones(m)
        col[0] = col[0] + m / 2
        rc = cp.outer(row, col)
        row = n / 2 * cp.ones(n)
        row[0] = row[0] + n / 2
        col = m * cp.cos(cp.pi * cp.arange(0, m) / m) - 2 * m
        col[0] = col[0] - m
        rc = rc + cp.outer(row, col)
        row = n / 2 * cp.ones(n)
        row[0] = n
        col = m / 2 *cp.ones(m)
        col[0] = m
        rc = cp.expand_dims(-1*rc + lamb*(h**2)*cp.outer(row, col),axis = 0).astype(cp.float64)
        y1 = cp.fft.fft(cp.concatenate([h*h*b, cp.zeros_like(b, dtype=cp.complex128)], axis=1),axis = 1)
        y1 = cp.real(y1[:,:m, :] * cp.expand_dims((cp.cos(cp.arange(0, m).reshape(-1, 1) * cp.pi / (2 * m)) -\
                                                 cp.sin(cp.arange(0, m).reshape(-1, 1) * cp.pi / (2 * m)) * 1j).reshape(-1,1),axis = 0))
        y1 = cp.fft.fft(cp.concatenate([y1.transpose(0,2,1).conj(), cp.zeros_like(y1.transpose(0,2,1), dtype=cp.complex128)], axis=1),axis = 1)
        y1 = cp.real(y1[:,:n, :] * cp.expand_dims((cp.cos(cp.arange(0, n).reshape(-1, 1) * cp.pi / (2 * n)) -\
                                                   cp.sin(cp.arange(0, n).reshape(-1, 1) * cp.pi / (2 * n)) * 1j).reshape(-1,1),axis = 0))
        # y1[:, 0, 0] = 0
        y1 = (y1 / rc).transpose(0,2,1).conj()
        # y1[: ,0, 0] = 0

        y2 = cp.real( cp.fft.fft(cp.concatenate([ cp.expand_dims(((cp.cos(cp.arange(0, m).reshape(-1, 1) * cp.pi / (2 * m)) ) -\
                                                   (cp.sin(cp.arange(0, m).reshape(-1, 1) * cp.pi / (2 * m))) * 1j),axis = 0)*y1,\
                                                 cp.zeros((bs, m, n))], axis=1),axis = 1))
        y2 = y2[:,:m,:].transpose(0,2,1).conj()
        y2 = cp.real( cp.fft.fft(cp.concatenate([cp.expand_dims((cp.cos(cp.arange(0, n).reshape(-1, 1) * cp.pi / (2 * n)) ),axis = 0)*y2 -\
                                                 cp.expand_dims((cp.sin(cp.arange(0, n).reshape(-1, 1) * cp.pi / (2 * n)) ),axis = 0)*y2*1j,\
                                                 cp.zeros((bs, n, m))], axis=1),axis = 1))
        return y2[:,:n, :].transpose(0,2,1).conj()
### Local-Intensity-Cluster的主函数    
def LIC_main(image,K,params):
    # image只要是[batch,width,height]的经过归一化的numpy数组即可
    # params配置好了参数
    image = torch.tensor(image,dtype = torch.float32).unsqueeze(1).cuda()
    phi = level_set_initialize(image)
#     print(phi.shape)
    for k in trange(params['maxepoch'],desc = 'Local-Intensity-Cluster'):
        with torch.no_grad(): 
            grad = cp.asarray(getGrad(phi,image,K,params).squeeze(1).detach().cpu().numpy())
            phi_new = cp.asarray(phi.squeeze(1).detach().cpu().numpy(),dtype = cp.float32)
#             print(grad.shape,grad.dtype)
            phi_new = fft_batch_helmholtz((phi_new - params['eta']*grad)/params['eta']/params['mu'],1.0/params['eta']/params['mu'])
            phi = torch.Tensor(phi_new.get()).unsqueeze(1).cuda()
    return phi
    
    
if __name__ == "__main__":
    
#     ### 使用example
#     X = np.random.rand(5,234,243)
#     X = torch.tensor(X,dtype = torch.float32).cuda()
#     K = create_gauss_kernel(1.0,5)
#     b = torch.ones_like(X).cuda()
#     image = torch.rand_like(b).cuda()
#     phi = LIC_main(image,K,params)
    
#     ### 检查初始mask:
#     phi = level_set_initialize(image)
#     plt.clf()
#     plt.contourf(phi[0,0,:,:].detach().cpu().numpy())
#     plt.show()

    ### 对于CT图像分割
    filename = "LIC_INP"
    imgL = []
    for idx in range(120,130): 
          fp = os.path.join(filename,"ser004img00{:d}.dcm".format(idx))
          image = np.expand_dims(rescale0_1(read_single_dicom_by_sitk(fp)),axis = 0) * 255.0
          imgL.append(image)
    imgL = np.squeeze(np.stack(imgL),axis = 1)
    # print(imgL.shape)
    # import pdb;pdb.set_trace()
    K = create_gauss_kernel(5.0,21)
    # print(image.shape)
    phi = LIC_main(torch.tensor(imgL,dtype = torch.float32),K,params)
    # import pdb;pdb.set_trace()
    plt.clf()
#     plt.contourf(phi[0,0,:,:].detach().cpu().numpy())
    for idx in range(120,130):
         plt.imshow(imgL[idx - 120,:,:])
         plt.contour(phi[idx - 120,0,:,:].detach().cpu().numpy(), levels=[0], colors='r')
         plt.savefig(os.path.join("LIC_RES","{:d}_.png".format(idx)))
         plt.clf()
         # plt.show()
    