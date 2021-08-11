%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          TimeNxtStep.m                              >
%  >                Set of Code that performs the executive role         >
%  >               of controlling the solutions at each time step        >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function [res,dt,itr] = TimeNxtStep(fid,Fluid,Rock,NumOutputs,Well,grid,Numerical,gridConne,...
                       Trans,Trans_Bd,res_Prev,dt,depth,steps,finT,dim,ph)
%Setting necessary intermediates used for calculations
N = grid.N;
CellVol = grid.dx*grid.dy*grid.dz;
res_Int = res_Prev;
%if(steps == 1)
%    res_Int.Vars = res_Prev.Vars.*1.000001;
%end
toler = Numerical.tol;
MaxIt = Numerical.maxIter;
t_Prev = res_Int.time;
x = 1;

%% SOLVERS 
%% Lagging Coefficient
            if strcmp(Numerical.solver,'lag') || strcmp(Numerical.solver,'Lag')
% Initializing the loop variables
                res_Int = ResUpdate(res_Int.Vars,Fluid,grid,Rock,dt,t_Prev,res_Int,0);
                itr = 1;
                dtredns = 0;  
                while(itr <= MaxIt && dtredns <= 10)
% Calculating Fluid Transmissibilities and adding relative 
                    [FTrans_I_W,FTrans_I_O] = FluidTrans(grid,Fluid,res_Int,gridConne,Trans,Trans_Bd);
        
% Combining Transmissibility due to Fluids and Gridding Structure
                     Trans_I_O    = FTrans_I_O; %.* Trans;
                     Trans_I_O_Bd = FTrans_I_O; %.* Trans_Bd;
                     Trans_I_W    = FTrans_I_W; %* Trans;
                     Trans_I_W_Bd = FTrans_I_W ;%.* Trans_Bd;
        
% Creating a Matrix comprising Transmissibilities based on Calculations        
                     Trans_Matrix_I = TransMatrix(grid,Trans_I_O,Trans_I_O_Bd,Trans_I_W,...
                           Trans_I_W_Bd,gridConne,N,ph);         
% Construct Accumulation Matrix         
                     [accum] = ConstructAccum( grid,Rock,Fluid,res_Int,res_Prev,ph )  ;
% Making the Vectors for Source/Sinks
        %[Q,res_Int] = SSConstruct(Well,res_Int,grid,Trans_I_Bd,ph);
                     [Q,res_Int] = SSConstruct(Well,res_Int,grid,0,ph,Trans_Matrix_I);
% Calculating the contribution of Gravity
%         gravity = 1/144 .* res_Int.Dens_O(:,1).*depth(:,1);
%         Grav    =  Trans_Matrix_I * gravity;
%         Grav = sparse(Grav);
        SolnsPrev = [preIterPressure, preIterSw];     
        SolnsPrev = reshape(SolnsPrev(:,1:2)',2*grid.N,1);
% Final Equation for material Balance
                     Left  = Trans_Matrix_I - accum./dt;
        %SLeft = sparse(Left); 
                    % Right = (-1)*accum./dt*res_Prev.Vars + Q; % + Grav;
        Right = (-1)*accum./dt*SolnsPrev + Q; % + Grav;
        %SRight = sparse(Right);
%% Method 1: Applying backslah (Gaussian Elimination)
                    Solution = Left\Right;
        %Solution = full(Solution);
% Method 2: BiCGStab Using ILU
%         tol = 1.0e-8; maxit = 1000;
%         setup.type = 'nofill';
%         [L1,U1]=ilu(sparse(Left),setup);
%         [Solution, flag, relres, itr] = bicgstab(sparse(Left), Right, tol ,maxit, L1, U1);
% Updating the intermediate solution used so far
                    res_Int = ResUpdate(Solution,Fluid,grid,Rock,dt,t_Prev,res_Int,1);
% Evaluating Convergence        
                    errorP = norm(res_Int.Press - preIterPressure)/norm(res_Int.Press);
                    errorS = norm(res_Int.Sat_W - preIterSw)/norm(res_Int.Sat_W);
                    if errorP < toler && errorS < toler
                        fprintf('>>> CONVERGED:  At TimeStep %i after %i Iterations with an Error of %d \n',steps,itr,errorP);
                        break;
                    elseif (errorP > toler ||errorS > toler) && itr >= MaxIt
                        dt = dt/2;
                        dtredns = dtredns + 1;
                        fprintf('!!!!! CONVERGENCE not achieved at TIMESTEP (%i) after Iteration Limit of (%i)......\n',steps,MaxIt);
                        fprintf('Time Step Reduction %i at ( %i  ,  %i )>>>>>>>>>>>>> Dt reduced to %d \n\n',dtredns,steps,itr,dt);
                        itr = 0;
                        x = 0;
                    end
                    if(NumOutputs.Itinf == 1 && x~=0)
                        fprintf('.... ITERATING:  At [%i, %i] , Error %d, %d \n',steps,itr,errorP, errorS);
                    end
                    itr = itr + 1;
                    preIterPressure = res_Int.Press;
                    preIterSw       = res_Int.Sat_W;
        
                end
                if (dtredns==0 && itr < Numerical.DtIncrCrit)||((finT - res_Int.time) < dt)
                    if((finT - res_Int.time) < dt)
                        dt = (finT - res_Int.time);
                        if(dt == 0)
                        return
                        end
                    else
                        if(~((finT - res_Int.time) < 2*dt))
                            dt = 2*dt;
                        end
                    end
                fprintf('!!!!! CONVERGENCE achieved at (%i,%i) without TS Reduction within %i Iterations \n',steps,itr,Numerical.DtIncrCrit);
                fprintf('>>> [ %i , %i ] Dt doubled to %d for next iteration............. at Time %d \n\n',steps,itr,dt,res_Int.time);
                end

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
 end

end


