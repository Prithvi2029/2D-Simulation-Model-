%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%  >     PETE 656 MATLAB FLOW SIMULATOR 1.0 - Final Project Spring 2021  >
%  >                          W/O Prithvi Singh Chauhan                  >
%  >                             09/05/2021                              >
%  >                            ResUpdate.m                              >
%  >                Set of Code that updates the PVT Properties          >
%  >                        at the incremented conditions                >
%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%
function res_Int = ResUpdate(Solution,Fluid,grid,Rock,dt,t_Prev,res_Int,p)
% Calculating the new properties at incremented pressure
if(p==1)
res_Int.time  = t_Prev + dt;
end
res_Int.Vars  = full(Solution);
res_Int.Press = Solution(1:2:end);
res_Int.Sat_W = Solution(2:2:end); 
%res_Int.Sat_W(res_Int.Sat_W > 0.7) = 0.7;
%res_Int.Sat_W(res_Int.Sat_W < 0.25) = 0.25;

res_Int.Sat_O = 1 - res_Int.Sat_W;
for k = 1:grid.Nz
    for j = 1:grid.Ny
        for i = 1:grid.Nx
            index = (k-1)*grid.Nx*grid.Ny + (j-1)*grid.Nx + i ;
% Incrementing the fluid properties and porous medium properties to the new
% conditions
            %Oil
%             res_Int.Visc_O = polyval(Fluid.viscosityFit,res_Int.Press);
%             res_Int.FVF_O = polyval(Fluid.FVFFit,res_Int.Press);
            res_Int.Visc_O(index,1)  = (Fluid.Oil.Visc_Cf(1)*...
                                        res_Int.Press(index) + Fluid.Oil.Visc_Cf(2));
            res_Int.FVF_O(index,1)   = (Fluid.Oil.FVF_Cf(1)*...
                                        res_Int.Press(index)  + Fluid.Oil.FVF_Cf(2));
            res_Int.Dens_O(index,1)  = (Fluid.Oil.RefDens_Oil*exp(Fluid.Oil.CompO*...
                                       (res_Int.Press(index) - Fluid.Oil.Ref_P)));
            res_Int.Poro(index,1)    = (Rock.poro(i,j,k)*exp(Rock.CompR*...
                                       (res_Int.Press(index) - Fluid.Oil.Ref_P)));
            if(res_Int.Sat_O(index,1) < 0.30) 
                res_Int.RelP_O(index,1) = 0.0 ;
            elseif (res_Int.Sat_O(index,1) > 0.75)
                res_Int.RelP_O(index,1) = 0.7 ;
            else
                res_Int.RelP_O(index,1) = 0.7 * ((res_Int.Sat_O(index,1)- 0.30)/(1 - 0.30 -0.25))^3;
            end
            %Water
            res_Int.Visc_W(index,1)  = (Fluid.H2O.Ref_Visc)*exp(Fluid.H2O.CompW*...
                                       (res_Int.Press(index) - Fluid.H2O.Ref_P));
            res_Int.Dens_W(index,1)  = (Fluid.H2O.RefDens_H2O*exp(Fluid.H2O.CompW*...
                                       (res_Int.Press(index) - Fluid.H2O.Ref_P)));
            res_Int.FVF_W(index,1)   = Fluid.H2O.RefFVF_H2O*exp(-Fluid.H2O.CompW*...
                                       (res_Int.Press(index) - Fluid.H2O.Ref_P));
            if(res_Int.Sat_W(index,1) < 0.25) 
                res_Int.RelP_W(index,1) = 0.0 ;
            elseif (res_Int.Sat_W(index,1) > 0.70)
                res_Int.RelP_W(index,1) = 0.08 ;
            else
                res_Int.RelP_W(index,1) = 0.08 * ((res_Int.Sat_W(index,1)- 0.25)/(1 - 0.30 -0.25))^2;
            end
        end
    end
%

end
end

  
                          

