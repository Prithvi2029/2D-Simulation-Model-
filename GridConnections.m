%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          GridConnections.m                          >
%  >                Code set to assign grid indices and                  >
%  >                  identify the nearest neighbours                    >
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function gridConne = GridConnections(grid);
%>>>>>>>>>>>>>>>>>>>>Determine maximum number of neighbours
if(grid.Ny == 1 && grid.Nz == 1)
    max_neighbours = 2;
elseif(grid.Nz == 1)
    max_neighbours = 4;
else
    max_neighbours = 6;
end
gridConne.MaxN = max_neighbours;
%% Cells are numbers along x followed by y and then along z directions.
for k = 1:grid.Nz
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            index = (k-1)*grid.Nx*grid.Ny + (j-1)*grid.Nx + i ;
            % Now at every index finding the number of neighbours
            for n = 1:max_neighbours
                % According to the numbering scheme
                % if  n = 1, West neighbour
                %     n = 2, East neighbour
                %     n = 3, South neighbour
                %     n = 4, North neighbour
                %     n = 5, Bottom neighbour
                %     n = 6, Top neighbour
                % These conditional statements discounts blocks on the boundary
                di = 0; dj = 0; dk = 0;
                if n == 1 && i == 1
                    continue;
                elseif n == 2 && i == grid.Nx
                    continue;
                elseif n == 3 && j == 1
                    continue;
                elseif n == 4 && j == grid.Ny
                    continue;
                elseif n == 5 && k == 1
                    continue;    
                elseif n == 6 && k == grid.Nz
                    continue;
                end
                % this loop identifies the index difference between the
                % gridblock and its neighbour
                if(n == 1)
                    di = -1;
                elseif(n == 2)
                    di = 1;
                elseif(n == 3)
                    dj = -1;
                elseif(n == 4)
                    dj = 1;
                elseif(n == 5)
                    dk = -1;
                elseif(n == 6)
                    dk = 1;
                end
                in = i + di;
                jn = j + dj;
                kn = k + dk;
                nindex = (kn - 1)*grid.Nx*grid.Ny + (jn - 1)*grid.Nx + in;
                % Assigning indexes of the neighbouring element for each
                % cell
                gridConne.neighbours(index,n)= nindex;
            end
            gridConne.i(index) = i;
            gridConne.j(index) = j;
            gridConne.k(index) = k;
        end
    end
end
