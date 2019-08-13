function [ EPat, W ] = PatEstimation( NL_mat, Self_arr, NSig, CurPat, Par )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code originated from 
% Weighted Nuclear Norm Minimization for Image Denoising, Version 1.0
% Shuhang Gu, Lei Zhang, Wangmeng Zuo, Xiangchu Feng
% https://github.com/csjunxu/WNNM_CVPR2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EPat = zeros(size(CurPat));
W    = zeros(size(CurPat));
for  i      =  1 : length(Self_arr)                                 % For each keypatch group
    Temp    =   CurPat(:, NL_mat(1:Par.patnum,i));                  % Non-local similar patches to the keypatch
    M_Temp  =   repmat(mean( Temp, 2 ),1,Par.patnum);
    Temp    =   Temp-M_Temp;
    
    E_Temp 	=   WNNM( Temp, Par.c, NSig, M_Temp, Par.reWeiIter); % WNNM Estimation
    EPat(:,NL_mat(1:Par.patnum,i))  = EPat(:,NL_mat(1:Par.patnum,i))+E_Temp;
    W(:,NL_mat(1:Par.patnum,i))     = W(:,NL_mat(1:Par.patnum,i))+ones(Par.patsize*Par.patsize,size(NL_mat(1:Par.patnum,i),1));
end
end
