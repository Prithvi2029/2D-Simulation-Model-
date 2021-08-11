%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          GroupedResults.m                           >
%  >                Code set to group results into a solution matrix     >
%  >                            at the every time step                   >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function group_res = GroupedResults(group_res,res,Well)

group_res.time = [group_res.time res.time];
group_res.elapset = [group_res.elapset res.elapset];
group_res.Press = [group_res.Press res.Press];
group_res.Visc_O = [group_res.Visc_O res.Visc_O];
group_res.FVF_O = [group_res.FVF_O res.FVF_O];
group_res.Poro = [group_res.Poro res.Poro];
group_res.Dens_O = [group_res.Dens_O res.Dens_O];
group_res.Visc_W = [group_res.Visc_W res.Visc_W];
group_res.FVF_W = [group_res.FVF_W res.FVF_W];
group_res.Dens_W = [group_res.Dens_W res.Dens_W];
group_res.Sat_W = [group_res.Sat_W res.Sat_W];
group_res.Sat_O = [group_res.Sat_O res.Sat_O];

% Source/Sink Rates
WellIdent = fieldnames(Well);
NumSS = length(WellIdent);
for i = 1:NumSS
    name = char(WellIdent(i));
    group_res.(name).oil_rate = [group_res.(name).oil_rate res.(name).oil_rate];
    group_res.(name).water_rate = [group_res.(name).water_rate res.(name).water_rate];
    group_res.(name).BHP = [group_res.(name).BHP res.(name).BHP];
end
end
    
    
    
    