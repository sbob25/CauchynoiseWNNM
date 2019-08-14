function [ energy ] = ComputeEnergy( ref_idx, ref_noisy_pat, rec_img, Par, lambda, gamma )

rec_pat = zeros(Par.patsize*Par.patsize,length(ref_idx));
k = 0;
for i  = 1:Par.patsize
    for j  = 1:Par.patsize
        k = k+1;
        E_patch = rec_img(i:end-Par.patsize+i,j:end-Par.patsize+j);
        E_patch = E_patch(:)';
        rec_pat(k,:) =  E_patch(:,ref_idx);
    end
end

[~,sing_mat,~] = svd(double(rec_pat),'econ');
energy = trace(sing_mat) + 0.5*lambda*sum(sum(log(gamma^2 + (rec_pat-ref_noisy_pat).^2)));

end
