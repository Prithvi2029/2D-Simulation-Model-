%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          NumJacobianComp.m                          >
%  >                Set of Code to calculate the Jacobi Matrix           >
%  >                      Numerically for the Residuals                  >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function JacobiM = NumJacobianComp(Trans_Matrix,Pressure_Prev,Pressure,accum,dt,...
                    Q,Grav,Numerical)
 
MatSz = size(Trans_Matrix,1);
prtbn = Numerical.PerturbJ;
JacobiM = [];
for i = 1:MatSz
                Pressure_Pt      = Pressure;
                Pressure_Pt(i,1) = Pressure_Pt(i,1) +  prtbn; 
                
% forming the sparse matricces saving memory and computational time
                Trans_Matrix   = sparse(Trans_Matrix);
                Pressure_Prev  = sparse(Pressure_Prev);
                Pressure_Pt    = sparse(Pressure_Pt);
                Pressure       = sparse(Pressure);
                accum          = sparse(accum);
                
% Calculating the residual Matrix
                ResidualM      = Trans_Matrix * Pressure - accum * ...
                                 (Pressure - Pressure_Prev)./dt - Q - Grav;
% Calculating the residual Matrix after Perturbations       
                PertdResidualM = Trans_Matrix * Pressure_Pt - accum * ...
                                 (Pressure_Pt - Pressure_Prev)./dt - Q - Grav;
                J = (PertdResidualM - ResidualM)./prtbn;
                JacobiM = [JacobiM J];
end
end            
