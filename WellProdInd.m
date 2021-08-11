%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          WellProdInd.m `````                        >
%  >                Code set to calculate the productivity index         >
%  >                              of a well                              >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function Well = WellProdInd(Rock,Well,grid)
%%Calculate the productivity index related to the well
% add it to the Well Object
WellIdent = fieldnames(Well);
NumSS = length(WellIdent);
for i = 1:NumSS
    %assuming well outer drainage radius is related to the grid size by a
    %relation of e^(-pi/2)
    name = char(WellIdent(i));
    ro = exp(-pi/2)*grid.dx;
    h = grid.dz;
    rw = Well.(name).rw;
    s  = Well.(name).skin;
    for j = 1:Well.(name).KInd
        wellperm(j) = Rock.perm(Well.(name).Wloci,Well.(name).Wlocj,j);
        Well.(name).WPI(j) = 0.001127*2*pi*wellperm(j)*h/(log(ro/rw) + s);
    end
end
end