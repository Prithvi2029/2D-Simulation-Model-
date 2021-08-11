%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          FluidTrans.m                               >
%  >                Set of Code that's Called in TimeNXTSTEP.m           >
%  >                for calculating fluid transmissibilities             >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function [FTrans_W,FTrans_O] = FluidTrans(grid,Fluid,res_Int,gridConne,Trans,Trans_Bd)

% Allocating maximum size to array

FTrans_W = zeros(grid.N,gridConne.MaxN);
FTrans_O = zeros(grid.N,gridConne.MaxN);
AA = grid.Agr;
BB = grid.Bgr;
for j=1:2:gridConne.MaxN
    if j==1
	    % West side blocks
		%OIL
		PWend = AA*res_Int.Press;
		DelP  = (res_Int.Press - PWend);
		DelP(1:grid.Ny:end) = 0;
		KrOW    = AA*res_Int.RelP_O;
		FVF_OW  = AA*res_Int.FVF_O ;
		Visc_OW = AA*res_Int.Visc_O ;
		FTOil1  = Trans(:,j).*res_Int.RelP_O./res_Int.FVF_O.*res_Int.Visc_O;		
		FTOil2  = Trans(:,j).*KrOW./FVF_OW.*Visc_OW;
		%WATER
		KrWW    = AA*res_Int.RelP_W;
		FVF_WW  = AA*res_Int.FVF_W ;
		Visc_WW = AA*res_Int.Visc_W ;
		FTWater1= Trans(:,j).*res_Int.RelP_W./res_Int.FVF_W.*res_Int.Visc_W;
		FTWater2= Trans(:,j).*KrWW./FVF_WW.*Visc_WW;
		% Upwinding Transmissibility according to flow
		for i = 1:length(DelP)
		    if DelP(i) >=0
				FTrans_O(i,j) = FTOil1(i);
				FTrans_W(i,j) = FTWater1(i);
			else
				FTrans_O(i,j) = FTOil2(i);
				FTrans_W(i,j) = FTWater2(i);
			end
		end
	elseif j==3
	    % North side blocks
		%OIL
		PNend = BB*res_Int.Press;
		DelP  = (res_Int.Press - PNend);
		DelP(1:grid.Nx:end) = 0;
		KrON    = BB*res_Int.RelP_O;
		FVF_ON  = BB*res_Int.FVF_O ;
		Visc_ON = BB*res_Int.Visc_O ;
		FTOil1  = Trans(:,j).*res_Int.RelP_O./res_Int.FVF_O.*res_Int.Visc_O;		
		FTOil2  = Trans(:,j).*KrON./FVF_ON.*Visc_ON;
		%WATER
		KrWN    = BB*res_Int.RelP_W;
		FVF_WN  = BB*res_Int.FVF_W ;
		Visc_WN = BB*res_Int.Visc_W ;
		FTWater1= Trans(:,j).*res_Int.RelP_W./res_Int.FVF_W.*res_Int.Visc_W;
		FTWater2= Trans(:,j).*KrWN./FVF_WN.*Visc_WN;
		% Upwinding Transmissibility according to flow
		for i = 1:length(DelP)
		    if DelP(i) >=0
				FTrans_O(i,j) = FTOil1(i);
				FTrans_W(i,j) = FTWater1(i);
			else
				FTrans_O(i,j) = FTOil2(i);
				FTrans_W(i,j) = FTWater2(i);
			end
		end
	end
end	