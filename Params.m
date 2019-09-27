function  [Par]=Params(nSig)
Par.nSig      = nSig;              % Variance of the noise image
Par.searchWin = 30;                % Non-local patch searching window
Par.c         = 23;                % Constant num for the weight vector
Par.innerloop = 2;                 % InnerLoop Num of between re-blockmatching
Par.reWeiIter = 3;
Par.patsize   = 6;                 % Patch size
Par.patnum    = 70;                % Initial Non-local Patch number
Par.iter      = 2;
Par.eta       = 0.6;               % Noise estimete parameter
Par.step      = 4;
Par.threshold = 1e-3;              % terminate condition of ADMM algorithm
