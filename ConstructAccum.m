%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                            Accum.m                                  >
%  >                Set of Code that performs the executive role         >
%  >               of controlling the solutions at each time step        >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function [ accum ] = ConstructAccum( grid,Rock,Fluid,res_Int,res_Prev,ph )
% Calculating Accumulation Terms 
         CellVol = grid.dx*grid.dy*grid.dz;
%          if(res_Int.Press == res_Prev.Press)
            Cr = Rock.CompR;
            Bw = Fluid.H2O.CompW;
            Bo = Fluid.Oil.CompO;
%          else
%              Cr = (res_Int.Poro-res_Prev.Poro)./(res_Int.Press - res_Prev.Press);
%              Bw = (1./res_Int.FVF_W - 1./res_Prev.FVF_W)./(res_Int.Press - res_Prev.Press);
%              Bo = (1./res_Int.FVF_O - 1./res_Prev.FVF_O)./(res_Int.Press - res_Prev.Press);
%          end
% Arranging Terms
        D11 = CellVol/5.6145.*((Cr./res_Int.FVF_W+Bw.*res_Int.Poro).*res_Int.Sat_W);
        D12 = CellVol/5.6145.*(res_Int.Poro./res_Int.FVF_W);
        D21 = CellVol/5.6145.*(res_Int.Sat_O.*(Cr./res_Int.FVF_O+res_Int.Poro.*Bo));
        D22 = CellVol/5.6145.*(-res_Int.Poro./res_Int.FVF_O);
        
        DiagA  = [D11, D22];
        DiagA  = reshape(DiagA(:,1:2)', 2*grid.N,1);
        
        ADiagL = [D21,zeros(grid.N,1)];
        ADiagL = reshape(ADiagL(:,1:2)',2*grid.N,1);

        ADiagU = [D12,zeros(grid.N,1)];
        ADiagU = reshape(ADiagU(:,1:2)',2*grid.N,1);
% Creating Sparse Matricces
        DiagASp  = sparse(1:ph*grid.N,1:ph*grid.N,DiagA(1:end),ph*grid.N,ph*grid.N);
        ADiagLSp = sparse(2:ph*grid.N,1:ph*grid.N-1,ADiagL(1:end-1),ph*grid.N,ph*grid.N);
        ADiagUSp = sparse(1:2*grid.N-1,2:ph*grid.N,ADiagU(1:end-1),ph*grid.N,ph*grid.N);
% Calculating Accumulation Term        
        accum = DiagASp + ADiagLSp + ADiagUSp;
%% Newton Raphson Method (Fully Implicit)
% elseif strcmp(Numerical.solver,'nr') || strcmp(Numerical.solver,'NR')
% % Setting necessary intermediates used for calculations
%     itr = 1;
%     dtredns = 0;
%     Pressure_Prev = res_Prev.Press;
% % Setting an initial guess in the pressure  
%     Pressure = Pressure_Prev-0.001;
% % Setting to run only till maximum Iterations or maximum number of time step reductions
%     while(itr <= MaxIt && dtredns <= 10)
% % Calculating Fluid Transmissibilities
%         FTrans_I = FluidTrans(grid,Fluid,res_Int,gridConne);
%         
% % Combining Transmissibility due to Fluids and Gridding Structure
%         Trans_I    = FTrans_I .* Trans;
%         Trans_I_Bd = FTrans_I .* Trans_Bd;
%         
% % Creating a Matrix comprising Transmissibilities based on Calculations        
%         Trans_Matrix_I = TransMatrix(grid,Trans_I,Trans_I_Bd,gridConne,N); 
% 
% % Calculating Accumulation Terms 
%         accum_m = (res_Int.Poro(:,1) .* Fluid.CompO ./ Fluid.RefFVF_Oil(:,1) ...
%                   + Rock.Ref_por * Rock.CompR ./ res_Prev.FVF_O(:,1))/5.6145;
%         accum   = diag(accum_m).*CellVol;
%         accum   = sparse(accum);
% 
% % Making the Vectors for Source/Sinks
%         [Q,res_Int] = SSConstruct(Well,res_Int,grid,Trans_I_Bd);
% 
% % Calculating the contribution of Gravity
%         gravity = 1/144 .* res_Int.Dens_O(:,1).*depth(:,1);
%         Grav    =  Trans_Matrix_I * gravity;
%         Grav = sparse(Grav);  
%         Trans_Matrix_I = sparse(Trans_Matrix_I);
%         
% % Forming the sparse matricces for Pressure
%         Pressure_Prev  = sparse(res_Prev.Press);
%         Pressure       = sparse(Pressure);
%         
% % Calculating the residual Matrix for the terms set up
%         ResidualM      = Trans_Matrix_I * Pressure - accum * ...
%                         (Pressure - Pressure_Prev)./dt - Q - Grav;
%         
% % Calculate the Jacobi Matrix for the terms set up
%         JacobiM = NumJacobianComp(Trans_Matrix_I,Pressure_Prev,Pressure,accum,dt,...
%                     Q,Grav,Numerical);
% % Incrementing the Pressures for next iteration
%         DelPress = (-1)*JacobiM\ResidualM;
%         Pressure = Pressure + full(DelPress);
% % Updating the intermediate solution used so far
%         res_Int = ResUpdate(Pressure,Fluid,grid,Rock,dt,t_Prev,res_Int);
% % Evaluating Convergence        
%         error = norm(Pressure-Pressure_Prev)/norm(Pressure);
%         if error < toler
%             fprintf(fid,'>>> CONVERGED:  At TimeStep %i after %i Iterations with an Error of %d \n',steps,itr,error);
%             break;
%         elseif error > toler && itr >= MaxIt
%             dt = dt/2;
%             dtredns = dtredns + 1;
%             fprintf(fid,'!!!!! CONVERGENCE not achieved at TIMESTEP (%i) after Iteration Limit of (%i)......\n',steps,itr);
%             fprintf(fid,'Time Step Reduction %i at ( %i  ,  %i )>>>>>>>>>>>>> Dt reduced to %d \n\n',dtredns,steps,itr,dt);
%             itr = 0;
%         end
%         if(NumOutputs.Itinf == 1 && x~=0)
%         fprintf(fid,'.... ITERATING:  At [%i, %i] , Error %d \n',steps,itr,error);
%         end
%         itr = itr + 1;
%         Pressure_Prev = full(Pressure);
%     end
%     res = res_Int;
%     if dtredns==0 && itr < Numerical.DtIncrCrit
%         if((finT - res_Int.time) < dt)
%             dt = (finT - res_Int.time);
%             if(dt == 0)
%                 return
%             end
%         else
%             if(~((finT - res_Int.time) < 2*dt))
%                 dt = 2*dt;
%             end
%         end
%             
%         fprintf(fid,'!!!!! CONVERGENCE achieved at (%i,%i) without TS Reduction within %i Iterations \n',steps,itr,Numerical.DtIncrCrit);
%         fprintf(fid,'>>> [ %i , %i ] Dt doubled to %d for next iteration............. at Time %d \n\n',steps,itr,dt,res_Int.time);
%     end
% else
%     fprintf('Invalid Solver Option, Available Options are:\n');
%     fprintf('1.Lagging Coefficient "Lag"\n');
%     fprintf('2.Newton Raphson      "NR"\n');
%         
% end



