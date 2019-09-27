function [ EPat, counting ] = WNNMDenoising( NL_mat, NSig, CurPat, Par)
patch_num = Par.patnum;
patch_size = Par.patsize;
reweight_iter = Par.reWeiIter;
C = Par.c;
patch_size2 = patch_size*patch_size;
ref_len = size(NL_mat, 2);

Y_stack = zeros(patch_size2, patch_num, ref_len);
meanY_stack = zeros(patch_size2, patch_num, ref_len);
X_stack = zeros(patch_size*patch_size, patch_num, ref_len);
EPat = zeros(size(CurPat));
counting = zeros(size(CurPat));

for i = 1:ref_len
    cur_stack = CurPat(:, NL_mat(1:patch_num,i));
    meanY_stack(:,:,i) = repmat(mean(cur_stack,2),1, patch_num);
    Y_stack(:,:,i) = cur_stack - meanY_stack(:,:,i);
end
parfor i = 1:ref_len
    [U, DY, V] = svd(Y_stack(:,:,i), 'econ');
    DDX = sqrt(max(diag(DY).^2 - patch_num*NSig^2, 0));
    for ii = 1:reweight_iter
        W_Vec = min(NSig^2, (C*sqrt(patch_num)*NSig^2)./(DDX+eps));
        DX = sign(DY).*max(abs(DY)-diag(W_Vec),0);
        if ii ~= reweight_iter
            DDX = diag(DX);
        else
            X_stack(:,:,i) = U*DX*V';
        end
    end
end
for i = 1:ref_len
    EPat(:,NL_mat(1:patch_num,i)) = EPat(:,NL_mat(1:patch_num,i)) + X_stack(:,:,i) + meanY_stack(:,:,i);
    counting(:,NL_mat(1:patch_num,i)) = counting(:,NL_mat(1:patch_num,i))+ones(patch_size*patch_size, patch_num);
end
end

