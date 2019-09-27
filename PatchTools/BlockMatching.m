function  [init_index] = BlockMatching(X, Par, neighbor_arr, num_arr, self_index_arr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code originated from 
% Weighted Nuclear Norm Minimization for Image Denoising, Version 1.0
% Shuhang Gu, Lei Zhang, Wangmeng Zuo, Xiangchu Feng
% https://github.com/csjunxu/WNNM_CVPR2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(num_arr);
init_index = zeros(Par.patnum,L);

parfor  i  =  1 : L
    patch = X(:,self_index_arr(i));
    neighbors = X(:,neighbor_arr(1:num_arr(i),i));
    dist = sum((repmat(patch,1,size(neighbors,2))-neighbors).^2);
    [~, index] = sort(dist);
    init_index(:,i)=neighbor_arr(index(1:Par.patnum),i);
end
