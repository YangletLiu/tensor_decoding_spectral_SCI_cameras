%TEST_GAPWNNM_TV_BIRD24 Test GAP-WNNM with GAP-TV for initialization.
%   See also GAPDENOISE_CACTI, DEMO_GAPDENOISE, GAPDENOISE, 
%     TEST_GAPDENOISE.
% clear; clc;
% [0] environment configuration
addpath(genpath('./algorithms')); % algorithms
addpath(genpath('./packages')); % packages
addpath(genpath('./utils')); % utilities
 
datasetdir = './dataset'; % dataset
resultdir  = './results'; % results

% [1] load dataset
para.type   = 'cassi'; % type of dataset, cassi or cacti
para.name   = 'toy31'; % name of dataset
para.number = 32; % number of frames in the dataset
datapath = 'toy31_cassi.mat';
% datapath = sprintf('%s/%s/%s%d_%s.mat',datasetdir,para.type,para.name,...
%     para.number,para.type);

if exist(datapath,'file')
    load(datapath); % mask, meas, orig (and para)
else
    error('File %s does not exist, please check dataset directory!',...
        datapath);
end

 meas = meas(161:300,161:300);
 mask = mask(161:300,161:300,:);
 orig = orig(161:300,161:300,:);


 %meas = meas(201:300,201:300);
 %mask = mask(201:300,201:300,:);
 %orig = orig(201:300,201:300,:);


para.nframe = 1; % number of coded frames in this test
para.MAXB   = 255;

[nrow,ncol,nmask] = size(mask);
%[nrow,ncol,nmask] = size(mask);
nframe = para.nframe; % number of coded frames in this test
MAXB = para.MAXB;

% [2] apply GAP-Denoise for reconstruction
% yall = meas/MAXB;

para.Mfunc  = @(z) A_xy(z,mask);
para.Mtfunc = @(z) At_xy_nonorm(z,mask);

para.Phisum = sum(mask.^2,3);
para.Phisum(para.Phisum==0) = 1;

% common parameters
para.lambda   =     1; % correction coefficiency
para.acc      =     1; % enable acceleration
para.flag_iqa =  true; % disable image quality assessments in iterations

%% [2.1] GAP-TV(-acc)
para.denoiser = 'tv'; % TV denoising
para.maxiter  = 500; % maximum iteration
para.tvweight = 5; % weight for TV denoising
para.tviter   = 5; % number of iteration for TV denoising
 
[vgaptv,psnr_gaptv,ssim_gaptv,tgaptv] = ...
     gapdenoise_cacti(mask,meas,orig,[],para);
 
 fprintf('GAP-%s mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
     upper(para.denoiser),mean(psnr_gaptv),mean(ssim_gaptv),tgaptv);
%% [2.6] GAP-WNNM-TV
para.denoiser = 'wnnm'; % WNNM denoising
  para.wnnm_int_fwise = true; % enable GAP-WNNM integrated (with frame-wise denoising)
    para.blockmatch_period = 20; % period of block matching
  para.sigma   = [12]/MAXB; % noise deviation (to be estimated and adapted)
  para.vrange  = 1; % value range
  para.maxiter = [8];
  para.patchsize = 32; % patch size
  para.iternum = 1; % iteration number in WNNM
  para.enparfor = 0; % enable parfor
  if para.enparfor % if parfor is enabled, start parpool in advance
      delete(gcp('nocreate')); % delete current parpool
      mycluster = parcluster('local');
      div = 1;
      while nmask/div > mycluster.NumWorkers
          div = div+1;
      end
      parnum = max(min(ceil(nmask/div),mycluster.NumWorkers),1); 
      parpool(mycluster,parnum);
  end
 %load abc.mat
 %vgaptv = vgapwnnm_tv;
load mask.mat
[vgapwnnm_tv,psnr_gapwnnm_tv,ssim_gapwnnm_tv,tgapwnnm2,psnrall_wnnmtv] = ...
    gapdenoise_cacti(mask,meas,orig,vgaptv,para);

tgapwnnm_tv = tgapwnnm2 + tgaptv;
delete(gcp('nocreate')); % delete parpool

fprintf('GAP-%s-TV mean PSNR %2.2f dB, mean SSIM %.4f, total time % 4.1f s.\n',...
    upper(para.denoiser),mean(psnr_gapwnnm_tv),mean(ssim_gapwnnm_tv),tgapwnnm_tv);

%% [3] save as the result .mat file and run the demonstration code
       %demo_gapdenoise.
 matdir = [resultdir '/savedmat'];
 if ~exist(matdir,'dir')
     mkdir(matdir);
 end
 save([matdir '/gapdenoise_1902293.mat']);
 save([matdir '/test_gapwnnm_tv_hsi' ...
     '_patchsize' num2str(para.patchsize) ...
     '_' para.name num2str(nframe*nmask) ...
     '_' char(datetime('today','format','yyMMdd')) ...
     '_' char(java.lang.System.getProperty('user.name')) '.mat']);
