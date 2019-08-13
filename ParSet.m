function  [Par]=ParSet(nSig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code modified from 
% Weighted Nuclear Norm Minimization for Image Denoising, Version 1.0
% Shuhang Gu, Lei Zhang, Wangmeng Zuo, Xiangchu Feng
% https://github.com/csjunxu/WNNM_CVPR2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Par.nSig      = nSig;              % Variance of the noise image
Par.searchWin = 30;                % Non-local patch searching window
Par.c         = 21;                % Constant num for the weight vector
Par.innerloop = 2;                 % InnerLoop Num of between re-blockmatching
Par.reWeiIter = 3;
Par.patsize   = 6;                 % Patch size
Par.patnum    = 70;                % Initial Non-local Patch number
Par.iter      = 2;
Par.eta       = 0.6;               % Noise estimete parameter
Par.step      = 2;
