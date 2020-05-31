function [E_Img PSNR SSIM]   =  NGmeet_Reconstruction1( y,opt)
 A  = @(x) opt.Mfunc(x);
 At = @(z) opt.Mtfunc(z); 
 Phisum = opt.Phisum;
 lambda = 0.1;
 y1 = zeros(size(y),'like',y);

rec_im0 = opt.v0;
rec_im           =    255*double(rec_im0);
AtY              =    At(y);
% beta             =    par.beta;  %0.01;

[h,h,w]           =    size(rec_im);  % size of the 2D image

Par   = ParSetH(5,31);
N_Img = reshape(rec_im,[h h w]);
delta =2;
E_Img            = N_Img;                                                         % Estimated Image
[Height, Width, Band]  = size(E_Img);  %% size of the 3D image
N = Height*Width;
TotalPatNum      = (Height-Par.patsize+1)*(Width-Par.patsize+1);         % Total Patch Number in the image
Average          = mean(N_Img,3);                      % Calculate the average band for fast spatial non-local searching
[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(Average, Par);   
% PreCompute all the patch index in the searching window 
for iter = 1 : Par.Iter
%First step: spectral dimension reduction 
   %% recontruction from single [0 1] data set
   E_Img = single(E_Img/255);
   yb = A(E_Img);
   y1 = y;
   E_Img = E_Img+lambda*(At((y1-yb)./Phisum)); % v=v+lambda*(At*A)^-1*At*dy
   E_Img = double(255*E_Img);
   %%     
   k_subspace = min(Par.k_subspace+delta*(iter-1),30);
   Y = reshape(E_Img, N, Band)';
   E_Img1 = E_Img;
   [E,~,~]= svd(Y,'econ');
    E=E(:,1:k_subspace);

    E_Img = reshape((E'*Y)', Height,Width, k_subspace);
    N_Img1 = reshape((E'*reshape(N_Img, N, Band)')', Height,Width, k_subspace); %%% add change N_Img1 as N_Img 
    Band1=k_subspace;
% %non-local patch grouping and noise estimation
    Average             =   mean(E_Img,3);
    [CurPat, Mat, Sigma_arr]	=	Cub2Patch( E_Img, N_Img1, Average, Par );

    if (mod(iter-1,2)==0)
        Par.patnum = max(Par.patnum - 10,round(Par.patnum/2));                                          % Lower Noise level, less NL patches
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

% reconstruction
%         b              =   AtY + beta * E_Img(:);
%        [X flag0]       =   pcg( @(x) Afun(x, At, A, beta), b, 0.5E-6, 400, [], [], E_Img(:));          
%         E_Img              =   reshape(X, Height, Width, Band); 
        
   
% time3 = toc
% estimation, can be ignored to speed up
    PSNR     =   csnr( reshape(E_Img,h*h,w),  255*reshape(opt.orig,h*h,w), 0, 0 );
%     fprintf( 'NGmeet CS Reconstruction, Iter %d : PSNR = %f \n', iter, PSNR );
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, NoiseLevel = %2.3f \n', iter, PSNR, sum(Sigma_arr)/TotalPatNum);
%     if iter<Par.Iter
%     E_Img = 0.1*N_Img+0.9*E_Img;
%     else
%     end
end
if isfield(opt,'orig')
    PSNR     =   csnr( reshape(E_Img,h*h,w), 255*reshape(opt.orig,h*h,w), 0, 0 );
    SSIM      =  cal_ssim( reshape(E_Img,h*h,w), 255*reshape(opt.orig,h*h,w), 0, 0 );
end
return;


%======================================================
% function  y  =  Afun(x, At, A, eta)
% y      =   At( A(x) ) + eta*x;  % eta * (Wei.*x);
% return;
% 
% 
% %====================================================================
% % Compressive Image Recovery Using DCT 
% %--------------------------------------------------------------------
% function  Rec_im0    =   DCT_CIR( y, par, A, At )
% ori_im      =    par.ori_im;
% [h w]       =    size(ori_im);
% im          =    At( y );
% im          =    reshape(im,[h w]);
% 
% lamada      =    1.5;  % 1.8, 1.2-1.7
% b           =    par.win*par.win;
% D           =    dctmtx(b);
% 
% for k   =  1:1
%     f      =   im;
%     for  iter = 1 : 300 %% 300   
%         
%         if (mod(iter, 50) == 0)
%             if isfield(par,'ori_im')
%                 PSNR     =   csnr( f, par.ori_im, 0, 0 );                
%                 fprintf( 'DCT Compressive Image Recovery, Iter %d : PSNR = %f\n', iter, PSNR );
%             end
%         end
%         
%         for ii = 1 : 3
%             fb        =   A( f(:) );
%             f         =   f + lamada.*reshape(At( y-fb ), h, w);
%         end        
%         f          =   DCT_thresholding( f, par, D );
%     end
%     im     =  f;
% end
% Rec_im0   =  im;
% return;

