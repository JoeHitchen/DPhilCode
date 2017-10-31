function [IntPos,Profiles] = DiffGrow_PDE_Routine(Stefan,LqdusGrad,LqdusInf,Biot,PDE,varargin)
% This function performs the PDE integration for a time-dependent run.
% It can be used to solve for the time-evolution of mushy sea ice, under a
% given set of conditions.
%
% Inputs:  Stefan - Stefan number
%          LqdusGrad - Liquidus gradient
%          LqdusInf  - Initial liquidus temperature
%          Biot      - Biot number
%          PDE       - PDE structure
%                      Required fields:
%                        zGrid         - Vertical grid
%                        tGrid         - Time grid
%                        MushyDynamics - Mushy (true) or simplified (false)
%                        RelTol        - Relative tolerance of iteration
%                        AbsTol        - Absolute tolerance of iteration
% (Opt.)   Ratios    - Thermal property ratio structure
%                      Possible fields:
%                        HeatCap   - Specific heat capacity ratio
%                        ThermCond - Thermal conductivity ratio
%
% Outputs: IntPos   - The time-dependent interface position
%          Profiles - The time-dependent temperature profiles
%
% (22/03/16)
    
    % Check for optional ratio inputs
    if nargin > 5
        
        % Set ratio structure
        Ratios = varargin{1};
    end
    
    % Set thermal conductivity ratio correction function
    if exist('Ratios','var') && isstruct(Ratios) && isfield(Ratios,'ThermCond') && Ratios.ThermCond ~= 1
        ThermCondCorrFcn = @(LqdFrac) (Ratios.ThermCond-1)*(1-LqdFrac);
    else
        ThermCondCorrFcn = @(LqdFrac) 0;
    end
    
    % Set heat capacity ratio correction function
    if exist('Ratios','var') && isstruct(Ratios) && isfield(Ratios,'HeatCap') && Ratios.HeatCap ~= 1
        HeatCapCorrFcn = @(LqdFrac) (Ratios.HeatCap-1)*(1-LqdFrac);
    else
        HeatCapCorrFcn = @(LqdFrac) 0;
    end
    
    % Set compositional Stefan number
    StefComp = Stefan/LqdusGrad;

    % Call PDE routine
    if PDE.MushyDynamics
        Profiles = pdepe(0, @MushyThermDiffEqn,@InitialState,@BoundaryConditions,PDE.zGrid,PDE.tGrid,odeset('RelTol',PDE.RelTol,'AbsTol',PDE.AbsTol));
    else
        Profiles = pdepe(0,@SimpleThermDiffEqn,@InitialState,@BoundaryConditions,PDE.zGrid,PDE.tGrid,odeset('RelTol',PDE.RelTol,'AbsTol',PDE.AbsTol));
    end
    Profiles = Profiles';
    
    % Calculate time-dependent interface positions
    IntPos = NaN(size(PDE.tGrid));
    for it = find(PDE.tGrid~=0)
        LqdInd       = find(Profiles(:,it)>LqdusInf); % Interpolating across the boundary makes no difference
        [~,tmpInd]   = unique(Profiles(LqdInd,it));
        LqdInd       = LqdInd(tmpInd);
        IntPos(it) = interp1(Profiles(LqdInd,it),PDE.zGrid(LqdInd),LqdusInf,'pchip');
    end
    IntPos(IntPos < 0) = NaN;
    
    
    %% PDE Subfunctions
    
    % Thermal diffusion equation function
    function [Mass,Flux,Source] = MushyThermDiffEqn(z,t,T,dTdz)            %#ok<INUSL>
        
        % Mushy region
        if T <= LqdusInf
            
            % Calculate liquid fraction
            LqdFrac = LqdusGrad./(LqdusGrad + LqdusInf - T);
            
            % Get mass and flux terms
            Mass = 1 + HeatCapCorrFcn(LqdFrac) + StefComp*LqdFrac.^2;
            Flux = (1+ThermCondCorrFcn(LqdFrac))*dTdz;
            
        % Liquid region
        else
            Mass = 1;
            Flux = dTdz;
        end
        Source = 0;
    end
    
    % Thermal diffusion equation function
    function [Mass,Flux,Source] = SimpleThermDiffEqn(z,t,T,dTdz)           %#ok<INUSL>
        if T <= LqdusInf
            Mass = 1 + StefComp;
        else
            Mass = 1;
        end
        Flux = dTdz;
        Source = 0;
    end
    
    % Initial condition function
    function T0 = InitialState(z)                                          %#ok<INUSD>
        T0 = 1;
    end
    
    % Boundary condition function
    function [pt,qt,pb,qb] = BoundaryConditions(zt,Tt,zb,Tb,t)             %#ok<INUSL,INUSD>
        pt = Tt;
        qt = -1/Biot;
        pb = Tb-1;
        qb = 0;
    end
end
