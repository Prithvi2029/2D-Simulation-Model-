%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          TransMatrix.m                              >
%  >               Set of Code that creates the Transmissibility         >
%  >                     matrix for the entire system                    >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function Trans_Matrix = TransMatrix(grid,Trans_G_O,Trans_GBd_O,Trans_G_W,...
                                     Trans_GBd_W,gridConne,N,ph)
Trans_Matrix = zeros(N*ph,N*ph);

switch gridConne.MaxN
	case 4
		%Transmissiblility of Oil 
		Trans_G_O(1:end-1,2) = Trans_G_O(2:end,1); %EW Boundary
		Trans_G_O(1:end-grid.Nx,4) = Trans_G_O(grid.Nx+1:end,3); %NS Boundary
		STmoil = -sum(Trans_G_O,2);
		
		%Transmissibility of Water
		Trans_G_W(1:end-1,2) = Trans_G_W(2:end,1); %EW Boundary
		Trans_G_W(1:end-grid.Nx,4) = Trans_G_W(grid.Nx+1:end,3); %NS Boundary
		STmH2O = -sum(Trans_G_W,2);		
		
		%Diagonals
		diag1     = [STmH2O,zeros(grid.N,1)];
		diag      = reshape(diag1(:,1:2)',ph*grid.N,1);
		Sparsedia = sparse(1:ph*grid.N,1:ph*grid.N,diag(1:end),ph*grid.N,ph*grid.N);
		
		L1dia     = [zeros(grid.N,1),STmoil];%non-generalized portion->zeros(grid.N,1) for Zero capillary pressures
		Ldia      = reshape(L1dia(:,1:2)',ph*grid.N,1);
		Ldia	  = sparse(2:ph*grid.N, 1:ph*grid.N-1, Ldia(2:end), ph*grid.N, ph*grid.N);
		
		U1dia     = [zeros(grid.N,1),Trans_G_O(:,2)];%non-generalized portion->zeros(grid.N,1) for Zero capillary pressures
		Udia	  = reshape(U1dia(:,1:2)',ph*grid.N,1);
		Udia	  = sparse(1:ph*grid.N-1,2:ph*grid.N,Udia(1:end-1),ph*grid.N,ph*grid.N);
		
		%West Flow Elements
		WfldiagW   = [Trans_G_W(:,1),zeros(grid.N,1)];
		WfldiagW   = reshape(WfldiagW(:,1:2)',ph*grid.N,1);
		WfldiagWSP = sparse(1+2:ph*grid.N,1:ph*grid.N-2,WfldiagW(1+2:end),ph*grid.N,ph*grid.N);
		WfldiagO   = [Trans_G_O(:,1),zeros(grid.N,1)];
		WfldiagO   = reshape(WfldiagO(:,1:2)',ph*grid.N,1);
		WfldiagOsp = sparse(1+3:ph*grid.N,1:ph*grid.N-3,WfldiagO(1+2:end-1),ph*grid.N,ph*grid.N);
		
		%East Band
		Efldiagsp  = sparse(1:ph*grid.N-2,1+2:ph*grid.N,WfldiagW(1+2:end),ph*grid.N,ph*grid.N);
		Efldiag	   = zeros(grid.N,2);
		Efldiag	   = reshape(Efldiag(:,1:2)',ph*grid.N,1);
		EfldiagUsp = sparse(1:ph*grid.N-3,1+3:ph*grid.N,Efldiag(1:end-3),ph*grid.N,ph*grid.N);
		
		%North Flow Elements
		NfldiagW   = [Trans_G_W(:,3),zeros(grid.N,1)];
		NfldiagW   = reshape(NfldiagW(:,1:2)',ph*grid.N,1);
		NfldiagWSp = sparse(1+ph*grid.Nx:ph*grid.N,1:ph*grid.N-ph*grid.Nx,NfldiagW(1+ph*grid.Nx:end),ph*grid.N,ph*grid.N);
		
		NfldiagO   = [Trans_G_O(:,3),zeros(grid.N,1)];
		NfldiagO   = reshape(NfldiagO(:,1:2)',ph*grid.N,1);
		NfldiagOSp = sparse(1+ph*grid.Nx+1:ph*grid.N,1:ph*grid.N-ph*grid.Nx-1,NfldiagO(1+ph*grid.Nx:end-1),ph*grid.N,ph*grid.N);
		
		NfldiagPc  = zeros(grid.N,2);
		NfldiagPc  = reshape(NfldiagPc(:,1:2)',ph*grid.N,1);
		NfldiagPcSp= sparse(1+ph*grid.Nx:ph*grid.N-1,2:ph*grid.N-ph*grid.Nx,NfldiagPc(1+ph*grid.Nx:end-1),ph*grid.N,ph*grid.N);
		
	   %South Flow
	   SfldiagSp   = sparse(1:ph*grid.N-ph*grid.Nx,1+ph*grid.Nx:ph*grid.N,NfldiagW(1+ph*grid.Nx:end),ph*grid.N,ph*grid.N);
	   SfldiagSpL  = sparse(2:ph*grid.N-ph*grid.Nx,1+ph*grid.Nx:ph*grid.N-1,NfldiagO(1+ph*grid.Nx:end-1),ph*grid.N,ph*grid.N);
	   SfldiagSpU  = sparse(1:ph*grid.N-ph*grid.Nx-1,2+ph*grid.Nx:ph*grid.N,NfldiagPc(1+ph*grid.Nx:end-1),ph*grid.N,ph*grid.N);
	   
	   Trans_Matrix = Sparsedia + Ldia + Udia + WfldiagWSP + WfldiagOsp + Efldiagsp + EfldiagUsp + NfldiagWSp ...
				      + NfldiagOSp + NfldiagPcSp + SfldiagSp + SfldiagSpL + SfldiagSpU;

    end
		
end