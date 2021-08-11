% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          GeometricTransm.m                          >
%  >                Code set to calculate the geometric                  >
%  >                  transmissibilty of a grid block                    >
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function [Trans,Trans_Bd] = GeometricTrans(gridConne,grid,Rock)
% Area along the faces of cells
AX = grid.dy * grid.dz;
AY = grid.dx * grid.dz;
AZ = grid.dx * grid.dy;

% Assigning sizes transmissibility matrices
Trans = zeros(grid.N, gridConne.MaxN);
Trans_Bd = zeros(grid.N, gridConne.MaxN);
% Assigning permeability for each grid block in vector form
for k = 1:grid.Nz
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            index = (k-1) * grid.Ny * grid.Nx + (j-1) * grid.Nx + i;
            perm(index) = Rock.perm(i,j,k);
        end
    end
end
% Calculating Transmissibility for each grid block
for k = 1:grid.Nz
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            index = (k-1) * grid.Ny * grid.Nx + (j-1) * grid.Nx + i;
            for m = 1:2:gridConne.MaxN
                nbor = gridConne.neighbours(index,m);
                if nbor ~= 0
                    if m == 1 || m == 2
                        Trans(index,m) = 2 * 0.001127 * AX * perm(index) * perm(nbor)/(perm(index)*grid.dx + perm(nbor)*grid.dx);
                    elseif m == 3 || m == 4
                        Trans(index,m) = 2 * 0.001127 * AY * perm(index) * perm(nbor)/(perm(index)*grid.dy + perm(nbor)*grid.dy);
                    elseif m == 5 || m == 6
                        Trans(index,m) = 2 * 0.001127 * AZ * perm(index) * perm(nbor)/(perm(index)*grid.dz + perm(nbor)*grid.dz);
                    end
                    
                elseif nbor == 0 && grid.PB(m) == 1
                    if m == 1 || m == 2
                        Trans_Bd(index,m) = 2 * 0.001127 * AX * perm(index)/grid.dx ;
                    elseif m == 3 || m == 4
                        Trans_Bd(index,m) = 2 * 0.001127 * AY * perm(index)/grid.dy;
                    elseif m == 5 || m == 6
                        Trans_Bd(index,m) = 2 * 0.001127 * AZ * perm(index)/grid.dz;
                    end                    
                end
            end
        end
    end
end
    
end

