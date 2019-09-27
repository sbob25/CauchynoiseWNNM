function [E_img]   =  XSubproblem( N_img, Par, NL_mat )

E_img = N_img;
[row, col] = size(N_img);

[neighbor_arr, num_arr, self_arr] =	NeighborIndex(N_img, Par);

for iter = 1 : Par.iter
    NSig = Par.eta*Par.nSig;
    CurPat = Im2Patch( E_img, Par );                      % image to patch and estimate local noise variance

    if (mod(iter-1,Par.innerloop)==0)
        tva = TotalVariation(E_img);
        tvb = TotalVariation(medfilt2(E_img));
        if tva/tvb < 3
            NL_mat = BlockMatching( CurPat, Par, neighbor_arr, num_arr, self_arr);% Caculate Non-local similar patches for each
        end
        NSig = Par.nSig;                       % First Iteration use the input noise parameter
    end

    [EPat, counting] = WNNMDenoising( NL_mat, NSig, CurPat, Par);
    E_img = Patch2Im( EPat, counting, Par.patsize, row, col );             
end
end


function [ v ] = TotalVariation(X)
    dx = [X(2:end,:) - X(1:end-1,:); zeros(1,size(X,2))];
    dy = [X(:,2:end) - X(:,1:end-1), zeros(size(X,1),1)];
    v = sum(sum(dx.^2+dy.^2));
end
