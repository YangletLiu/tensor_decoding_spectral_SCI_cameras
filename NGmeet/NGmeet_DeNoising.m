function [E_Img]= NGmeet_DeNoising( N_Img, O_Img, Par,iter)
% reference paper: Non-local Meets Global: An Integrated Paradigm for Hyperspectral Denoising
% oriData3_noise  N_Img        Input noisy 3-D image
% OriData3        O_Img        Reference image, for the PSNR computing of each step
% k_subspace  Par.k_subspace   The initial spectral rank, can be estimated by HySime
% delta=2                       iteration of k_subspace
% for pavia city and CAVE k_subspace = 5+2*(iter-1);
delta =1;
E_Img            = N_Img;                                                         % Estimated Image
[Height, Width, Band]  = size(E_Img);  
N = Height*Width;
TotalPatNum      = (Height-Par.patsize+1)*(Width-Par.patsize+1);         % Total Patch Number in the image
Average          = mean(N_Img,3);                      % Calculate the average band for fast spatial non-local searching
[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(Average, Par);   
% PreCompute all the patch index in the searching window 

% for iter = 1 : Par.Iter 
%First step: spectral dimension reduction 
   k_subspace = Par.k_subspace+delta*(iter-1);
   Y = reshape(E_Img, N, Band)';
   E_Img1 = E_Img;
%    [w Rw] = estNoise(Y,'additive');
%    Rw_ori = Rw;
%    Y = sqrt(inv(Rw_ori))*Y;
%    img_ori = reshape(Y', Height, Width, Band);
%    [w Rw] = estNoise(Y,'additive');
%    [~, E]=hysime(Y,w,Rw);
    [E,~,~]= ntsvd(Y,'econ');
%   [E,~,~]= svd(Y,'econ');
   
   
   
   
   
   
    E=E(:,1:k_subspace);

    E_Img = reshape((E'*Y)', Height,Width, k_subspace);
    N_Img1 = reshape((E'*reshape(N_Img, N, Band)')', Height,Width, k_subspace); %%% add change N_Img1 as N_Img 
    Band1=k_subspace;

% %non-local patch grouping and noise estimation
    Average             =   mean(E_Img,3);
    [CurPat, Mat, Sigma_arr]	=	Cub2Patch( E_Img, N_Img1, Average, Par );

    if iter>=1%(mod(iter-1,2)==0)
        Par.patnum = Par.patnum - 5;                                          % Lower Noise level, less NL patches
        NL_mat  =  Block_matching(Mat, Par, Neighbor_arr, Num_arr, Self_arr);  % Caculate Non-local similar patches for each
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr))*sqrt(k_subspace/Band);                      % First Iteration use the input noise parameter
        end
    end
%     time2=toc
% non-local low-rank denoising
    [Spa_EPat, Spa_W]    =  NLPatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par); 
% reconstruct patches to 3-D image
    [Spa_Img, Spa_Wei]   =  Patch2Cub( Spa_EPat, Spa_W, Par.patsize, Height, Width, Band1 );       % Patch to Cubic
    E_Img = Spa_Img./Spa_Wei;
    E_Img = reshape(reshape(E_Img, Height*Width, k_subspace)*E',Height,Width, Band);

% time3 = toc
% estimation, can be ignored to speed up
    [PSNR,SSIM,~,~] = evaluate(O_Img/255,E_Img/255,Height,Width);PSNR = mean(PSNR);SSIM = mean(SSIM);
%    PSNR = mean(0);SSIM = mean(0);
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, SSIM = %2.3f, NoiseLevel = %2.3f \n', iter, PSNR, SSIM, sum(Sigma_arr)/TotalPatNum);
    if iter<Par.Iter
    E_Img = 0.1*N_Img+0.9*E_Img;
    else
%     end
end
end





function [U,S,V] = ntsvd(A,fftOP,parOP)

% determine size of tensor
sa = size(A);
la = length(sa);
[n1,n2,n3]=size(A);
fl=0;


if ~exist('parOP','var')
    parOP = false;
end

if ~exist('fftOP','var')
   fftOP = false; 
end

%% Perform FFT along 3 to P axeses of Tensor


% conjugate symetric trick.
if la == 3
    if n2 > n1
        transflag=1;
        A=tran(A);
        nn1=n1;
        n1=n2;
        n2=nn1;
    end
    U = zeros(n1,n1,n3);
    S = zeros(n1,n2,n3);
    V = zeros(n2,n2,n3);
    % for P orders
else
    sU =sa;         %determines proper size
    sU(2) = sU(1);
    sV = sa;        %determines proper size
    sV(1) = sV(2);
    U = zeros(sU);  %pre allocated
    S = zeros(sa);
    V = zeros(sV);
end

for i = 3:la
    A = fft(A,[],i);
end

faces = prod(sa(3:la));     %determine # of faces

if la == 3
    % Do the conjugate symetric trick here.
    if isinteger(n3/2)
        endValue = int16(n3/2 + 1);
        [U, S, V] = takeSVDs(U,S,V,A,endValue,parOP);
        
        for j =n3:-1:endValue+1
            U(:,:,j) = conj(U(:,:,n3-j+2));
            V(:,:,j) = conj(V(:,:,n3-j+2));
            S(:,:,j) = S(:,:,n3-j+2);
        end

    else % if isinteger(n3/2)
        endValue = int16(n3/2 + 1);        
        [U,S,V] = takeSVDs(U,S,V,A,endValue,parOP);
        
        for j =n3:-1:endValue+1
            U(:,:,j) = conj(U(:,:,n3-j+2));
            V(:,:,j) = conj(V(:,:,n3-j+2));
            S(:,:,j) = S(:,:,n3-j+2);
        end
    end %if isinteger(n3/2)
%% for 4+ dimensional tensors do not perform the
% the conjugate trick.
else % if la == 3  
    [U, S, V] = takeSVDs(U,S,V,A,faces,parOP);
end

%%

if ~fftOP
    [U S V] = ifft_T(U,S,V);
end

if exist('transflag','var')
    Uold =U; U=V; S=tran(S); V=Uold;  
end


end
%% BEGIN SUBFUNCTIONS
%
%%
function [U, S, V] = takeSVDs(U,S,V,A,endI,runPar)

if ~exist('runPar','var')
    runPar = false;
end
    
if ~runPar || matlabpool('size') == 0

    for i=1:endI
        [U1,S1,V1]=svd(A(:,:,i));
        U(:,:,i)=U1; S(:,:,i)=S1; V(:,:,i)=V1;
    end
else
    
    parfor i=1:endI
        [U1,S1,V1]=svd(A(:,:,i));
        U(:,:,i)=U1; S(:,:,i)=S1; V(:,:,i)=V1;
    end
end

end

