%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                          ReadInputs.m                               >
%  >                Code set to read inputs in for Setting up            >
%  >                            Simulation                               >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function[dim,ph,Fluid, Rock, Well, grid, NumOutputs, Numerical,Media, depth] ...
                                                    = ReadInputs(Filename)
%
NumOutputs.OutT = 0;
%% Command lines to read the data from input file
fileID = fopen(Filename, 'r');
if fileID < 0
    disp('There was an error opening the file, Make sure it exists');
    %stop;
else
% Begin Reading input

    tline = fgetl(fileID);
    while (~strcmp(tline,'END'))
        C = textscan(tline,'%s');
        C = C{1};
% >>>>>>>>>>Read Problem Dimensions
        if strcmp(C{1},'Dimension')
           dim = str2double(C{2});
        elseif strcmp(C{1},'Phases')
           ph = str2double(C{2}); 
        elseif strcmp(C{1},'Depth')
           gendepth = C{2};
           if (~strcmp(gendepth,'gener'))
                P = C(2:length(C));
                depth = str2double(P);
                %Simple Entry of depths
           end
% >>>>>>>>>>Read Rock Properties
        elseif strcmp(C{1},'Rock_Compr') 
            Rock.CompR   = str2double(C{2});
        elseif strcmp(C{1},'Ref_Poro')
            Rock.Ref_por = str2double(C{2});
% >>>>>>>>>>Read Media Properties
        elseif strcmp(C{1},'RelPermEq')
            Media.RelP_eq = str2double(C{2});
        elseif strcmp(C{1},'Initial_P')
            Media.Int_P = str2double(C{2});
        elseif strcmp(C{1},'Oil_Saturation')
            Media.SatO = str2double(C{2});
            Media.SatW = 1 - Media.SatO;
% >>>>>>>>>>Read Fluid Properties
% >> Oil Properties
        elseif strcmp(C{1},'Oil')
            t=-1;
            while (t~=1)
                tline = fgetl(fileID);
                C = textscan(tline,'%s');
                C = C{1};
                if strcmp(C{1},'Reference_P')
                    Fluid.Oil.Ref_P = str2double(C{2});
                elseif strcmp(C{1},'Oil_Ref_Visc')
                    Fluid.Oil.Ref_Visc = str2double(C{2}); 
                 elseif strcmp(C{1},'Oil_Visc_Coeff')
                    Fluid.Oil.Visc_Cf = [str2double(C{2}) str2double(C{3})];
                elseif strcmp(C{1},'Oil_Compr')
                    Fluid.Oil.CompO = str2double(C{2});
                elseif strcmp(C{1},'Oil_FVF')
                    Fluid.Oil.RefFVF_Oil = str2double(C{2}); 
                elseif strcmp(C{1},'Oil_FVF_Coeff')
                    Fluid.Oil.FVF_Cf = [str2double(C{2}) str2double(C{3})];
                elseif strcmp(C{1},'Oil_Density')
                    Fluid.Oil.RefDens_Oil = str2double(C{2}); 
                else
                    t = 1;
                end
            end
% >> Water Properties
        elseif strcmp(C{1},'Water')
            t=-1;
            while (t~=1)
                tline = fgetl(fileID);
                C = textscan(tline,'%s');
                C = C{1};
                if strcmp(C{1},'Reference_P')
                    Fluid.H2O.Ref_P = str2double(C{2});
                elseif strcmp(C{1},'H2O_Ref_Visc')
                    Fluid.H2O.Ref_Visc = str2double(C{2}); 
                elseif strcmp(C{1},'H2O_Compr')
                    Fluid.H2O.CompW = str2double(C{2});
                elseif strcmp(C{1},'H2O_FVF')
                    Fluid.H2O.RefFVF_H2O = str2double(C{2}); 
                elseif strcmp(C{1},'H2O_Density')
                    Fluid.H2O.RefDens_H2O = str2double(C{2}); 
                else
                    t = 1;
                end
            end
% >>>>>>>>>>Read Well Properties
        elseif strcmp(C{1},'WellName')
            WellName = C{2};
            t = -1; 
            while (t~=1)
                tline = fgetl(fileID);
                C = textscan(tline,'%s');
                C = C{1};
                if strcmp(C{1},'rwell')
                    Well.(WellName).rw = str2double(C{2});
                elseif strcmp(C{1},'Prodtype')
                    Well.(WellName).control = C{2};
                elseif strcmp(C{1},'Welltype')
                    Well.(WellName).wt = C{2};
                elseif strcmp(C{1},'BHP')
                    Well.(WellName).BHP = str2double(C{2});
                elseif strcmp(C{1},'Rate')
                    Well.(WellName).RATE = str2double(C{2});
                elseif strcmp(C{1},'LocationI')
                    Well.(WellName).Wloci = str2double(C{2});
                elseif strcmp(C{1},'LocationJ')
                    Well.(WellName).Wlocj = str2double(C{2});
                elseif strcmp(C{1},'KIndex')
                    Well.(WellName).KInd = str2double(C{2});
                elseif strcmp(C{1},'skin')
                    Well.(WellName).skin = str2double(C{2});
                else
                    t = 1;
                end
            end
% >>>>>>>>>>Read Grid Properties
        elseif strcmp(C{1},'dx')
                grid.dx = str2double(C{2});
        elseif strcmp(C{1},'dy')
                grid.dy = str2double(C{2});
        elseif strcmp(C{1},'dz')
                grid.dz = str2double(C{2});
% >>>>>>>>>>Read Boundary Conditions
        elseif strcmp(C{1},'PB')
                grid.PB = [str2num(C{2});str2num(C{3});str2num(C{4});...
                            str2num(C{5});str2num(C{6});str2num(C{7})];
        elseif strcmp(C{1},'PBP')
                grid.PBP = [str2num(C{2});str2num(C{3});str2num(C{4});...        
                            str2num(C{5});str2num(C{6});str2num(C{7})];
% >>>>>>>>>>Read Numerical Conditions
        elseif strcmp(C{1},'tstep_int')
                Numerical.intdt = str2double(C{2});
        elseif strcmp(C{1},'t_final')
                Numerical.fint = str2double(C{2});
        elseif strcmp(C{1},'tol')
                Numerical.tol = str2double(C{2});
        elseif strcmp(C{1},'PerturbJ')
                Numerical.PerturbJ = str2double(C{2});
        elseif strcmp(C{1},'MaxIter')
                Numerical.maxIter = str2double(C{2});
        elseif strcmp(C{1},'solver')
                Numerical.solver = C{2};
        elseif strcmp(C{1},'DtIncrementCrit')
                Numerical.DtIncrCrit = str2double(C{2});
% >>>>>>>>>>Read Output Properties               
        elseif strcmp(C{1},'PrintOutputs')
                NumOutputs.PO = str2double(C{2});
        elseif strcmp(C{1},'PrintIteration')
                NumOutputs.Itinf = str2double(C{2});
        elseif strcmp(C{1},'PrintTimes')
                Ts = C(2:length(C));
                NumOutputs.OutT = str2double(Ts);
                %Simple Entry of Output Times
        end
        tline = fgetl(fileID);
    end
%%  Working on Input Files Read, Reading PoroPerm Data   
    if dim == 1
        data = load('data1D.mat');
        Rock.perm = data.perm.';
        Rock.poro = data.poro.';  
    elseif dim == 2
        data = load('data2D.mat');
        Rock.perm = data.perm';
        Rock.poro = data.poro;
    elseif dim == 3
        data = load('data.mat');
        Rock.perm = data.perm;
        Rock.poro = data.poro;
    end
%% Deciding the Number of Grid Blocks
    grid.Nx = size(Rock.perm,1);
    grid.Ny = size(Rock.perm,2);
    if dim==3
        grid.Nz = size(Rock.perm,3);
    else
        grid.Nz = 1;
    end
    grid.N = grid.Nx*grid.Ny*grid.Nz;
%% Deciding the Depth of System Blocks
    if strcmp(gendepth,'gener')
        depth = zeros(grid.N,1);
        j=1;
        for k=1:grid.Nz
            depth(j:k*grid.Nx*grid.Ny) = grid.dz*(grid.Nz-k)+grid.dz/2;
            j=j+grid.Nx*grid.Ny;
        end
     end
%% Boundary Conditions to be set
     if dim==1
        grid.PB=grid.PB(1:2,1);
        grid.PBP=grid.PBP(1:2,1);
     elseif dim==2
        grid.PB=grid.PB(1:4,1);
        grid.PBP=grid.PBP(1:4,1);
     elseif dim==3
        grid.PB=grid.PB(1:6,1);
        grid.PBP=grid.PBP(1:6,1);
     end
end
%Oil Properties
% Fluid.pRefo=[400,800,1200,1600,2000,2400,2800,3200,3600,4000,4400,...
%     4800,5200,5600];                                                       % Oil Reference Pressure
% Fluid.mu_o=[1.17,1.14,1.11,1.08,1.06,1.03,1,.98,.95,.94,.92,...
%     .91,.9,.89];                                                           % @Oil Ref. pressure - oil viscosity [cp]
% Fluid.B_o=[1.012,1.009,1.005,1.001,.996,.99,.988,.985,.98,.975,...
%     .97,.965,.96,.955];                                                    % @Oil Ref. pressure - oil formation volume factor [rb/stb]
Numerical.Schedule=[Numerical.intdt:0.1:10,10.2,10.6,11,11.5,12:1:279,280:0.1:284.9,285:0.5:290,291:1:366];
Numerical.tt=1;
Numerical.Time=Numerical.Schedule(Numerical.tt);
fclose(fileID);
end





















        