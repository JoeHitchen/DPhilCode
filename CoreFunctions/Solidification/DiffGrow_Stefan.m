function [InterfacePosition,OutputProfiles,Depths] = DiffGrow_Stefan(Stefan,LqdusInf,Comp,varargin)
% This function solves for the interface position of the Stefan problem for
% a perfectly conducting boundary and without salt diffusion by solving the
% implicit equation directly.
%
% Inputs:  Stefan    - Thermal Stefan number
%          LqdusInf  - Liquidus temperature at infinity
%
%          Comp  - Computation structure
%                    Required fields:
%                      DomainDepth - Computational domain depth
%                      NGrid       - Number of grid points in each domain
%                                    (Two values)
%
% (Opt.)   Ratio - Parameter ratios (Default values are 1)
%                    Possible fields:
%                      Dens - Density ratio
%                      Cond - Conductivity ratio
%                      Diff - Diffusivity ratio
%
% (Opt.)   Time  - Time for non-self similar solutions
%
% Outputs: InterfacePosition - Scaled position of the mush/liquid
%                                interface.
%          OutputProfiles    - Enthalpy, temperature, bulk salinity, liquid
%                                salinity and liquid fraction profiles.
%          Depths            - Depths associated with the output profiles
%
% (16/03/15)

    % Unpack ratios
    if nargin > 4 && isstruct(varargin{1})
        error('Input processing not written!')
    else
        Ratio.Dens = 1;
        Ratio.Cond = 1;
        Ratio.Diff = 1;
    end

    % Unpack time scaling
    if nargin > 5
        Time = varargin{2};
    else
        Time = 1;
    end
    
    % Set interface condition function
    G_Fcn = @(Var) sqrt(pi)*Var.*erf(Var).*exp(Var.^2);
    F_Fcn = @(Var) sqrt(pi)*Var.*erfcx(Var);
    IntCondFcn = @(Pos) Ratio.Cond/Ratio.Diff*LqdusInf./G_Fcn(Pos/(2*sqrt(Ratio.Diff))) - ...
        (1-LqdusInf)./F_Fcn(Pos/2) - Ratio.Dens*Stefan;
    
    % Find interface
    InterfacePosition = fzero(IntCondFcn,0.15);
    
    % Perform self-similarity inversion
    InterfacePosition = InterfacePosition*sqrt(Time);
    
    % Set domain depths
    Solid_Depths  = ChebGrid(Comp.NGrid(1),InterfacePosition);
    Liquid_Depths = InterfacePosition + ChebGrid(Comp.NGrid(2),Comp.DomainDepth-InterfacePosition);
    Depths  = [Solid_Depths;Liquid_Depths];
    
    % Set profiles
    OutputProfiles = NaN(length(Depths),5);
    OutputProfiles(1:Comp.NGrid(1)+1,2)   = LqdusInf*BiotCoolingProfile(Solid_Depths,Time,inf,Ratio.Diff)/BiotCoolingProfile(InterfacePosition,Time,inf,Ratio.Diff);
    OutputProfiles(Comp.NGrid(1)+2:end,2) = 1-(1-LqdusInf)*BiotCoolingProfile(Liquid_Depths,Time,inf,1)/BiotCoolingProfile(InterfacePosition,Time,inf,1);
    OutputProfiles(:,3) = ones(length(Depths),1);
    OutputProfiles(:,4) = [ones(size(Solid_Depths));zeros(size(Liquid_Depths))];
    OutputProfiles(:,5) = NaN(size(Depths));
    OutputProfiles(:,1) = OutputProfiles(:,2) - Stefan*(1-OutputProfiles(:,5));
    
    % Avoid duplication
    Depths         = Depths        ([1:Comp.NGrid(1)+1,Comp.NGrid(1)+3:end]);
    OutputProfiles = OutputProfiles([1:Comp.NGrid(1)+1,Comp.NGrid(1)+3:end],:);
end

