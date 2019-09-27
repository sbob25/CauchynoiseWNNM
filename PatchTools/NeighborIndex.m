function  [ Neighbor_arr, Num_arr, SelfIndex_arr ]  =  NeighborIndex(img, Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code modified from 
% Weighted Nuclear Norm Minimization for Image Denoising, Version 1.0
% Shuhang Gu, Lei Zhang, Wangmeng Zuo, Xiangchu Feng
% https://github.com/csjunxu/WNNM_CVPR2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function Precompute the all the patch indexes in the Searching window
% -Neighbor_arr is the array of neighbor patch indexes for each keypatch
% -Num_arr is array of the effective neighbor patch numbers for each keypatch
% -SelfIndex_arr is the index of keypatches in the total patch index array 
SW      	=   Par.searchWin;
s           =   Par.step;
TempR       =   size(img,1)-Par.patsize+1;
TempC       =   size(img,2)-Par.patsize+1;
R_GridIdx	=   (1:s:TempR);
R_GridIdx	=   [R_GridIdx R_GridIdx(end)+1:TempR];
C_GridIdx	=   (1:s:TempC);
C_GridIdx	=   [C_GridIdx C_GridIdx(end)+1:TempC];

Idx         =   (1:TempR*TempC);
Idx         =   reshape(Idx, TempR, TempC);
R_GridH     =   length(R_GridIdx);    
C_GridW     =   length(C_GridIdx); 

Neighbor_arr    =   int32(zeros((2*SW+1)*(2*SW+1), R_GridH*C_GridW));
Num_arr         =   int32(zeros(1,R_GridH*C_GridW));
SelfIndex_arr   =   int32(zeros(1,R_GridH*C_GridW));

for  i  =  1 : R_GridH
    for  j  =  1 : C_GridW    
        OffsetR     =   R_GridIdx(i);
        OffsetC     =   C_GridIdx(j);
        Offset1  	=  (OffsetC-1)*TempR + OffsetR;
        Offset2   	=  (j-1)*R_GridH + i;
                
        top         =   max( OffsetR-SW, 1 );
        button      =   min( OffsetR+SW, TempR );        
        left        =   max( OffsetC-SW, 1 );
        right       =   min( OffsetC+SW, TempC );     
        
        NL_Idx     =   Idx(top:button, left:right);
        NL_Idx     =   NL_Idx(:);

        Num_arr(Offset2)  =  length(NL_Idx);
        Neighbor_arr(1:Num_arr(Offset2),Offset2)  =  NL_Idx;   
        SelfIndex_arr(Offset2) = Offset1;
    end
end
