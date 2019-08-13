function  [X] =  WNNM( Y, C, nSig, m, iter )
    try
        [U,SigmaY,V] =   svd(full(Y),'econ');
    catch
        [U,SigmaY,V] =   svdecon(full(Y));
    end

    patch_num       = size(Y,2);
    Temp         =   sqrt(max( diag(SigmaY).^2 - patch_num*nSig^2, 0 ));

    for i=1:iter
        W_Vec    =   min(nSig^2, (C*sqrt(patch_num)*nSig^2)./( Temp + eps ));
        SigmaX   =  soft(SigmaY, diag(W_Vec));
       	Temp     = diag(SigmaX);
    end
    
    X =  U*SigmaX*V' + m;
return;

function y = soft(x,tau)

y = sign(x).*max(abs(x)-tau,0);