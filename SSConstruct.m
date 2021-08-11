%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          SSConstruct.m                              >
%  >                  Set of Code that creates the Vector                >
%  >                  representing the Source/Sink Terms                 >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function [Q,res] = SSConstruct(Well,res,grid,Trans_Bd,ph,Trans_Matrix_I)
% Initializing
WellIdent = fieldnames(Well);
NumSS = length(WellIdent);
Q = zeros(grid.N*ph,1);
PBSumQ = zeros(grid.N,1);

%% Running loop for each Source/Sink
for i = 1:NumSS
    name    = char(WellIdent(i));
    RateSumO = 0;
    RateSumW = 0;
    for j = 1:Well.(name).KInd
        index = (j-1)*grid.Nx*grid.Ny + (Well.(name).Wlocj-1)*grid.Nx ...
                 + Well.(name).Wloci;
        if strcmp(Well.(name).wt,'Prd')
% Constant Pressure Production             
            if strcmp(Well.(name).control,'Pressure')
                if Well.(name).BHP < res.Press(index)
                    % Oil Production
                    Q(2*index) = Well.(name).WPI(j) * (res.Press(index) ... 
                                                  - Well.(name).BHP);
%                     Well.(name).WellIndex=Well.(name).WPI*res.RelP_O(index)/(res.Visc_O(index)...
%                                                     *res.FVF_O(index));
%                     Q(2*index) = Well.(name).WPI(j) * ( ... 
%                                                   - Well.(name).BHP);
                    Q(2*index) = res.RelP_O(index)* Q(2*index)/(res.Visc_O(index)...
                                                    *res.FVF_O(index));
                    RateSumO   = RateSumO + Q(2*index);
                    %Trans_Matrix_I(2*index,2*index-1)= Trans_Matrix_I(2*index,2*index-1)-Well.(name).WellIndex;
                    % Water Production
                    Q(2*index-1) = Well.(name).WPI(j) * (res.Press(index) ... 
                                                  - Well.(name).BHP);
%                     Q(2*index-1) = Well.(name).WPI(j) * ( ... 
%                                                   - Well.(name).BHP);
                    Q(2*index-1) = res.RelP_W(index)* Q(2*index-1)/(res.Visc_O(index)...
                                                    *res.FVF_O(index));  
%                     Well.(name).WellIndex=Well.(name).WPI*res.RelP_W(index)/(res.Visc_W(index)...
%                                                     *res.FVF_W(index));
                    %Trans_Matrix_I(2*index,2*index-1)= Trans_Matrix_I(2*index-1,2*index-1)-Well.(name).WellIndex;                            
                    RateSumW   = RateSumW + Q(2*index-1);       
                    end
                BHP = Well.(name).BHP;
% Constant Rate Production            
            elseif strcmp(Well.(name).control,'Rate')
                    Q(2*index)   = Well.(name).RATE / Well.(name).KInd * ...
                                   res.RelP_O(index)/(res.RelP_O(index) + res.RelP_W(index));
                    Q(2*index-1) = Well.(name).RATE / Well.(name).KInd * ...
                                   res.RelP_W(index)/(res.RelP_O(index) + res.RelP_W(index));
                    BHP_O    =  Q(2*index)*res.Visc_O(index)*res.FVF_O(index)...
                                /Well.(name).WPI(j);
                    BHO_W    =  Q(2*index-1)*res.Visc_W(index)*res.FVF_W(index)...
                                /Well.(name).WPI(j);
                    BHP      = res.Press(index) - (BHP_O*Q(2*index)/(Q(2*index)+Q(2*index-1))+...
                                                 BHP_W*Q(2*index-1)/(Q(2*index)+Q(2*index-1)))
                    % Impossible case with -ve Bottom hole Pressure
                    if BHP < 0
                        BHP      = 0;
                        Q(2*index) = 0;
                        Q(2*index-1) = 0;
                    end
                    RateSumO  = RateSumO + Q(2*index);
                    RateSumW  = RateSumW + Q(2*index-1);
            end
        % Injection Well
        else
            % Constant Pressure Injection             
            if strcmp(Well.(name).control,'Pressure')
                if Well.(name).BHP < res.Press(index)
                    Q(2*index-1) = Well.(name).WPI(j) * (res.Press(index) ... 
                                                  - Well.(name).BHP);
                    Q(2*index-1) = res.RelP_W(index)* Q(index)/(res.Visc_O(index)...
                                                    *res.FVF_O(index));                    
                    RateSumW   = RateSumW + Q(2*index-1);
                end
                BHP = Well.(name).BHP;
            % Constant Rate Injection            
            elseif strcmp(Well.(name).control,'Rate')
                    Q(2*index-1) = -Well.(name).RATE / Well.(name).KInd;
                    BHP      = res.Press(index) - ...
                            Q(2*index-1)*res.Visc_W(index)*res.FVF_W(index)...
                            /Well.(name).WPI(j);
                    % Impossible case with -ve Bottom hole Pressure
                    if BHP < 0
                        BHP      = 0;
                        Q(2*index-1) = 0;
                    end
                    RateSumW  = RateSumW + Q(2*index-1);     
            end
        end
    end
    
    res.(name).oil_rate = RateSumO;
    res.(name).water_rate = RateSumW;
    res.(name).BHP    = BHP; 
end

% The Effect due to Pressure Boundary needs to be adjusted
% for k = 1:length(grid.PB)
%     PBTempQ = Trans_Bd(:,k).* grid.PBP(k);
%     PBSumQ =  PBSumQ + PBTempQ;
% end

%Q = Q - PBSumQ;
Q = sparse(Q);
end