%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                               Main.m                                >
%  >                    The MAIN FUNCTION is the branch                  >
%  >           that controls the running of all other functions          > 
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clc; clear; close all;
fid = fopen( 'result.txt', 'wt' );
%% Input Data
[dim,ph, Fluid, Rock, Well, grid, NumOutputs, Numerical, Media, depth] = ...
ReadInputs('InputPr5.txt');
N = grid.N;
toler = Numerical.tol;
MaxIt = Numerical.maxIter;
%% Initializing Geometric Attributes, i.e. Connectivity etc
gridConne = GridConnections(grid);
%% Accounting the contribution of the well, i.e. source/sink
% This is needed for cases where well is being produced at a constant
% pressure
Well = WellProdInd(Rock,Well,grid);
%% Calculating Geometric Transmissibilities
[Trans,Trans_Bd] = GeometricTrans(gridConne,grid,Rock);

%% Following Initialization of Source and Sink Properties and Block Transmissibilities
% Running Simulation Time Loop
% Initializing the input
[res,Fluid]  = InitialT(Fluid,Media,grid, Rock, Well); 

res_Prev  = res;

res_Prev.Vars = [res_Prev.Press, res_Prev.Sat_W];
res_Prev.Vars = reshape(res_Prev.Vars(:,1:2)',2*N,1);
group_res = res_Prev;
res_Int = res_Prev;
preIterPressure = res_Prev.Press; 
preIterSw       = res_Prev.Sat_W;
preIterVars     = res_Prev.Vars;
dt     = Numerical.intdt;
time   = res.time + dt;
%res_Int.Vars = res_Prev.Vars.*1.000001;

%%Creating the requisite matrices
grid.Agr = [diag(ones(grid.N,1)),zeros(grid.N,1)];
grid.Agr = grid.Agr(:,2:end);
grid.Agr(1:15:end,:) = 0;
grid.Agr = sparse(grid.Agr);
grid.Bgr = [zeros(15,grid.N);diag(ones(grid.N,1))];
grid.Bgr=sparse(grid.Bgr(1:end-15,:));

if NumOutputs.OutT(1)~=0
    Numerical.fint = [NumOutputs.OutT' Numerical.fint];
end

steps  = 1;
for i = 1:length(Numerical.fint)      
       while time < Numerical.fint(i)
            finT =  Numerical.fint(i);
%Setting necessary intermediates used for calculations
            t_Prev = res_Prev.time;
            x = 1;

%% SOLVERS 
%% Lagging Coefficient
            if strcmp(Numerical.solver,'lag') || strcmp(Numerical.solver,'Lag')
% Initializing the loop variables
                %res_Int = ResUpdate(res_Int.Vars,Fluid,grid,Rock,dt,t_Prev,res_Int,0);
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
%% Preconditioning using ILU
                    tol = 1.0e-4; maxit = 1000;
                    setup.type = 'nofill';
                   [L1,U1]=ilu(sparse(Left),setup);
%% Method 1: Applying backslash (Gaussian Elimination)
%                     Solution = Left\Right;
%% Method 2: BiCGStab Using preconditioner
                  [Solution, flag, relres, itr] = bicgstab(sparse(Left), Right, tol ,maxit, L1, U1);
%%           BiCGSTAB with no preconditioner
           %         [Solution, flag, relres, itr] = bicgstab(sparse(Left), Right, tol ,maxit);        
%% Method 3: CG with no preconditioner
%                    [Solution, flag, relres, itr] = pcg(sparse(Left), Right, tol, maxit);
%%           CG with preconditioner
%                    [Solution, flag, relres, itr] = pcg(sparse(Left), Right, tol, maxit, L1, U1);
%% Method 4: GMRES with no preconditioner
              %      [Solution, flag, relres, itr] = gmres(sparse(Left), Right, [], tol, maxit);
%%           GMRES with preconditioner
           %         [Solution, flag, relres, itr] = gmres(sparse(Left), Right, [], tol, maxit, L1, U1);
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
                        if(dt ==0)
                        break;
                        end
                    else
                        if(~((finT - res_Int.time) < 2*dt))
                            dt = 2*dt;
                    end
                end
                fprintf('!!!!! CONVERGENCE achieved at (%i,%i) without TS Reduction within %i Iterations \n',steps,itr,Numerical.DtIncrCrit);
                fprintf('>>> [ %i , %i ] Dt doubled to %d for next iteration............. at Time %d \n\n',steps,itr,dt,res_Int.time);
                end
            end   
            res = res_Int;
            group_res = GroupedResults(group_res,res,Well);
            time      = res.time;
            fprintf('>>> Printing DATA After ( %i  ,  %i ) ....................  TIME is %d days \n\n',steps,itr,time);
            res_Prev  = res;
            steps = steps + 1;
            %res_Int.Vars = res_Int.Vars.*1.000001;
        end
        dt = Numerical.intdt;
end
% k = 1;
% times = [0,30,60,90,120,150,180,210,240,270,300,330,360,365];
% for i = 1:(size(group_res.time,2))
%     if(k<=length(times))
%         if abs(group_res.time(i)- times(k))<=0.9 ||abs(group_res.time(i)- 30)<=0.01
%             Print1.time(k)  = group_res.time(i); 
%             Print1.Press(:,k) = group_res.Press(:,i);       
%             Print1.Sat_O(:,k)= group_res.Visc_O(:,i);
%             Print1.Sat_W(:,k)= group_res.Sat_W(:,i);
%             Print1.FVF_O(:,k) = group_res.FVF_O(:,i);
%             Print1.Poro(:,k)  = group_res.Poro(:,i);
%             WellIdent = fieldnames(Well);
%             NumSS = length(WellIdent);
%             for s = 1:NumSS
%                 name = char(WellIdent(s));
%                 Print1.(name).water_rate(k) = group_res.(name).water_rate(i); 
%                 Print1.(name).oil_rate(k) = group_res.(name).oil_rate(i); 
%                 Print1.(name).BHP(k)    = group_res.(name).BHP(i);  
%             end
%             k = k + 1;
%           end
%     end
% end
% if(dim == 2)
% X_Co_ord = 0:grid.dx:grid.dx*grid.Nx-grid.dx;
% Y_Co_ord = 0:grid.dy:grid.dy*grid.Ny-grid.dy;
% %figure(1)
%    % [XGrid, YGrid] = meshgrid(X_Co_ord,Y_Co_ord);
%    % pcolor(XGrid,YGrid,reshape(Print.Press(:,1),15,15))
%    % str = sprintf('Saturation distribution [psia] at Time %d days', Print.time(1));pbaspect([1 1 1]);
% 
% for t = 1:length(Print1.time) 
%     figure(2);
%     subplot(5,3,t)
%     [XGrid, YGrid] = meshgrid(X_Co_ord,Y_Co_ord);
%     contourf(XGrid,YGrid,reshape(Print1.Press(:,t),15,15),'ShowText','on')
%     str = sprintf('Pressure distribution [psia] at Time %d days', times(t));
%     title(str)
%     colormap(jet);h = colorbar;
%     hold on
%     ylabel(h, 'Pressure psia');caxis;pbaspect([1 1 1]);
% end
% hold off
% for l = 1:length(Print1.time) 
%     figure(3);
%     subplot(5,3,l)
%     [XGrid, YGrid] = meshgrid(X_Co_ord,Y_Co_ord);
%     pcolor(XGrid,YGrid,reshape(Print1.Sat_W(:,t),15,15))
%     str = sprintf('Saturation distribution [psia] at Time %d days', times(l));
%     title(str)
%     colormap(jet);h = colorbar;
%     hold on 
%     ylabel(h, 'Water Saturation');caxis;pbaspect([1 1 1])
% end
% hold off
% end
% 
% yy = smooth(group_res.time,group_res.('www1').water_rate,0.1,'loess');
% yy1=smooth(group_res.time,group_res.('www1').oil_rate,0.1,'loess');
% yyaxis left;
% figure(5)
% plot(group_res.time,yy);
% xlabel('time days');
% ylabel('Water rate');
% yyaxis right;
% plot(group_res.time,yy1);
% ylabel('Oil Rate');
% backax1 = gca;
% co1 = ax.Color;ax.Color = 'black';ax.LineWidth = 2;