# tensor_decoding_spectral_SCI_cameras
We apply tensor decoding scheme to decode SCI cameras.

This repository contains the code for the paper Exploiting Low-rank Tensor and Nonlocal Similarity for Decoding Hyperspectral SCI Cameras

Snapshot compressive imaging (SCI) cameras compress multiple spectral frames into a single measurement frame. An accurate and high-efficient decoding method for SCI cameras is important. However, existing algorithms such as GAP-TV and DeSCI are not satisfactory due to low accuracy and extensive executing time respectively. In this paper, we propose a new algorithm Low-rank Tensor and Nonlocal Similarity (LTNS) that addresses these two problems. 

Our proposed LTNS algorithm achieves an average improvement of 3.98dB PSNR over DeSCI and 5.63dB PSNR than GAP-TV. Also, our algorithm (LTNS) reduces the running time to less than one minute, comparing with the DeSCI’s 8 hours’ running time. We compare our algorithm with GMM-TP and MMLE-GMM as well.

This code uses the bird dataset and toy dataset, which are coded aperture snapshot spectral imaging (CASSI) dataset. The method we use to compare in the code is our LTNS and GAP-TV with several iterations.

## Usage
0. Requirements are MATLAB(R) with Parallel Computing Toolbox (parfor for multi-CPU acceleration).
1. Download this repository via git
```
git clone https://github.com/hust512/tensor_decoding_spectral_SCI_cameras
```
or download the [zip file](https://github.com/hust512/tensor_decoding_spectral_SCI_cameras/master.zip) manually.

2. Test the LTNS algorithm via
```matlab
Demo_bird_cassi_desci.m
```


## Structure of directories

| directory  | discription  |
| :--------: | :----------- | 
| `algorithms` | MATLAB function of the algorithms proposed in the article | 
| `packages`   | algorithms adapted from the state-of-art algorithms (adapted)|
| `dataset`    | Bird and toy (CASSI) dataset |
| `results`   | Save the result of reconstruction|
| `utils`      | utility functions |

## Platform
The test platform is MATLAB 2019b operating on macOS, also you can run on any machine with MATLAB(R) and Parallel Computing Toolbox, operating on Windows(R), Linux, or Mac OS. No GPU is needed to run this code.

Overall, our algorithm’s running time is around one minute and other algorithms except DeSCI range from ten seconds to forty minutes. DeSCI’s running time will be longer than eight hours.

## Reference
[1] W. He, Q. Yao, C. Li, N. Yokoya and Q. Zhao, “Non-local meets global: an integrated paradigm for hyperspectral denoising.” IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR), Long Beach, CA, USA, 2019, pp. 6861-6870.

[2] Y. Liu, X. Yuan, J. Suo, D. J. Brady and Q. Dai, “Rank minimization for snapshot compressive imaging,” in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 41, no. 12, pp. 2990-3006, 1 Dec. 2019.
https://github.com/hust512/DeSCI


## Contact
[X.-Y.Liu, Columbia University](http://www.tensorlet.com/ ) 

[Huping Ding, Columbia University](mailto:hd2436@columbia.edu ) 

[Tingyi Wang, Columbia University](mailto:tw2677@columbia.edu ) 
