function [ NL_mat, ref_idx, ref_noisy_pat ] = InitialBlockMatching(X, Y, Par)
XPat = Im2Patch(X, Par);
[neighbor_arr, num_arr, self_arr] =	NeighborIndex(X, Par);
NL_mat  =  BlockMatching( XPat, Par, neighbor_arr, num_arr, self_arr);

% To reduce computational cost, we compute energy only for the
% representative patch. We take the patch with the smallest variance of
% pixel values as the represntative patch.
fPat = XPat(:,NL_mat(1,:));
std_fPat = std(fPat,1);
[~,idx] = min(std_fPat);
ref_idx = NL_mat(:,idx);

ref_noisy_pat = zeros(Par.patsize*Par.patsize,length(ref_idx));
k = 0;
for i  = 1:Par.patsize
    for j  = 1:Par.patsize
        k = k+1;
        N_patch = Y(i:end-Par.patsize+i,j:end-Par.patsize+j);
        N_patch = N_patch(:)';
        ref_noisy_pat(k,:) = N_patch(:,ref_idx);
    end
end
end
