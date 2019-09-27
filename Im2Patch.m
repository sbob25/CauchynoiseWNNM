function  [ Pat_img ]  =  Im2Patch( E_img, Par )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code modified from 
% Weighted Nuclear Norm Minimization for Image Denoising, Version 1.0
% Shuhang Gu, Lei Zhang, Wangmeng Zuo, Xiangchu Feng
% https://github.com/csjunxu/WNNM_CVPR2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_pat_num = (size(E_img,1)-Par.patsize+1)*(size(E_img,2)-Par.patsize+1); 
Pat_img           =   zeros(Par.patsize*Par.patsize, total_pat_num);
k           =   0;
for i  = 1:Par.patsize
    for j  = 1:Par.patsize
              k     =  k+1;
        E_patch     =  E_img(i:end-Par.patsize+i,j:end-Par.patsize+j);
        Pat_img(k,:)      =  E_patch(:)';
    end
end
