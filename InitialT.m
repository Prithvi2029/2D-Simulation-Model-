%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          InitialT.m                                 >
%  >                Code set to initialize the time varying              >
%  >                     properties during simulation                    >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [res,Fluid] = InitialT(Fluid,Media, grid, Rock, Well)

%initial time

res.time = 0;
res.elapset = 0;

%% Initializing fluid properties
% Oil
res.Press   = ones(grid.N,1)* Media.Int_P;
res.Visc_O  = ones(grid.N,1)*(Fluid.Oil.Visc_Cf(1)*Media.Int_P + Fluid.Oil.Visc_Cf(2))  ;
res.FVF_O   = ones(grid.N,1)*(Fluid.Oil.FVF_Cf(1)*Media.Int_P + Fluid.Oil.FVF_Cf(2));
% % 2. initial oil viscosity
% Fluid.viscosityFit = polyfit(Fluid.pRefo,Fluid.mu_o,3);
% res.Visc_O = polyval(Fluid.viscosityFit,res.Press);
% % 3. initial oil FVF
% Fluid.FVFFit = polyfit(Fluid.pRefo,Fluid.B_o,3);
% res.FVF_O = polyval(Fluid.FVFFit,res.Press);
res.Dens_O  = ones(grid.N,1)*(Fluid.Oil.RefDens_Oil*exp(Fluid.Oil.CompO*(Media.Int_P - Fluid.Oil.Ref_P)));
res.Sat_O   = ones(grid.N,1)* Media.SatO;
if(Media.SatO <= 0.30) 
    res.RelP_O = ones(grid.N,1)*0.0 ;
elseif (Media.SatO >= 0.75)
    res.RelP_O = ones(grid.N,1)*0.7 ;
else
    res.RelP_O = ones(grid.N,1)* 0.7 * ((Media.SatO - 0.30)/(1 - 0.30 -0.25))^3;
end
% Water
res.Press   = ones(grid.N,1)* Media.Int_P;
res.Visc_W  = ones(grid.N,1)*(Fluid.H2O.Ref_Visc)  ;
res.Dens_W  = ones(grid.N,1)*(Fluid.H2O.RefDens_H2O*exp(Fluid.H2O.CompW*(Media.Int_P - Fluid.H2O.Ref_P)));
res.FVF_W   = ones(grid.N,1)*(Fluid.H2O.RefFVF_H2O*exp(-Fluid.H2O.CompW*(Media.Int_P - Fluid.H2O.Ref_P)));
res.Sat_W   = ones(grid.N,1)*(1 - Media.SatO);
% Calculating Relative Permeability of Water                   
if((1 - Media.SatO) <= 0.25)
    res.RelP_W = ones(grid.N,1)*0.0 ;
elseif ((1 - Media.SatO) >= 0.70 )
    res.RelP_W =ones(grid.N,1)* 0.08 ;
else
    res.RelP_W = ones(grid.N,1)* 0.08 * (((1 - Media.SatO) - 0.25)/(1 - 0.30 -0.25))^2;
end 
%initializing porous medium properties

for k = 1:grid.Nz
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            index    = (k-1) *grid.Nx*grid.Ny + (j-1) *grid.Nx + i;
            res.Poro(index,1) = Rock.poro(i,j,k)*exp(Rock.CompR*(Media.Int_P - Fluid.Oil.Ref_P));%Fluid.Ref_P));  
        end
    end
end

%Initializing Sources/Sink Properties

WellIdent = fieldnames(Well);
NumSS = length(WellIdent);
for i = 1:NumSS
    name = char(WellIdent(i));
    res.(name).prod_O = Well.(name).RATE;
    res.(name).BHP = Well.(name).BHP;
    res.(name).oil_rate = zeros(1,1);
    res.(name).water_rate = zeros(1,1);
end
end