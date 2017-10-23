function [InterfacePosition,OutputProfiles,OutputDepths,SearchVarFinal,MinOutput] = DiffGrow_Mushy_Cond_NS(Stefan,LqdusGrad,LqdusInf,Comp,varargin)
% This function solves for the interface position in the mushy model for a 
% perfectly conducting boundary and without salt diffusion using the 
% shooting method.
%
% Inputs: Stefan        - The Stefan number.
%         LqdusGrad     - The liquidus gradient.
%         LqdusInf      - The initial liquidus temperature
%         Comp          - Computation structure
%           Required fields:
%             DomainDepth - Computational domain depth
%             NGrid       - Number of grid points in each domain
%                            (Two values)
% (Opt.)  Initial guess - The initial guess of the three search variables.
%                           1. The interface position.
%                           2. The surface thermal gradient.
%                           3. The surface liquid fraction.
%
% Outputs: InterfacePosition - The interface position/growth rate.
%          OutputProfiles    - The enthalpy, thermal, bulk salinity, liquid
%                              salinity and liquid fraction profiles.
%          OutputDepths      - The depths corresponding to the calculated
%                              profiles.
%          SearchVarFinal    - The final values of the three search
%                              variables. For large Stefan numbers, it is
%                              often necessary to use the output from a
%                              smaller Stefan number as input.
%          MinOutput         - The value of the minimisation function.
%
% (14/01/15)
    
    % Check for initial guess
    if nargin > 4
        SearchVarInitial = varargin{1};
    else
        IntPosFcn = @(Pos) LqdusInf - erf(Pos/2);
        SearchVarInitial(1) = fzero(IntPosFcn,1);
        SearchVarInitial(2) = 1/sqrt(pi);
        SearchVarInitial(3) = 1/(1+LqdusInf/LqdusGrad);
    end
    
    % Unpack parameters
    DomainDepth = Comp.DomainDepth;
    
    % Create integration domain function
    IntDomFcn = @(SerVar) {[0,SerVar(1)];[SerVar(1),DomainDepth]};
    
    % Create initial condition function
    InitCondFcn = @(SerVar) [0,SerVar(2),SerVar(3)];
    
    % Create interface condtions functions
    InterfaceCondTemp = @(Profs) Profs(end,1:2);
    
    % Create residual condition function
    ResidualCond = @(Profs,IntInd) [Profs(end,1)-1,Profs(IntInd(1),1)-LqdusInf,Profs(IntInd(1),3)-1];
    
    % Perform minimisation
    options = optimoptions('fsolve','Display','off','TolFun',1e-6);
    SearchVarFinal    = fsolve(@MinimisationFunction,SearchVarInitial,options);
    InterfacePosition = SearchVarFinal(1);
    
    % Generate final profiles
    [Depths,ProfilesF,IntIdx] = GetProfiles(IntDomFcn(SearchVarFinal),InitCondFcn(SearchVarFinal));
    Depths    = Depths(   [1:IntIdx(1),(IntIdx(2)+1):end]);
    ProfilesF = ProfilesF([1:IntIdx(1),(IntIdx(2)+1):end],:);
        
    % Generate grid
    Mushy_Depths  = ChebGrid(Comp.NGrid(1),InterfacePosition);
    Liquid_Depths = InterfacePosition + ChebGrid(Comp.NGrid(2),Comp.DomainDepth-InterfacePosition);
    OutputDepths  = [Mushy_Depths;Liquid_Depths(2:end)];
    
    % Extrapolate profiles
    OutputProfiles(:,2) = interp1(Depths,ProfilesF(:,1),OutputDepths,'spline');
    OutputProfiles(:,3) = ones(Comp.NGrid*[1;1]+1,1);
    OutputProfiles(:,4) = max((LqdusInf-OutputProfiles(:,2))/LqdusGrad,0);
    OutputProfiles(:,5) = OutputProfiles(:,3)./(1+OutputProfiles(:,4));
    OutputProfiles(:,1) = OutputProfiles(:,2)-Stefan*(1-OutputProfiles(:,5));
    
    OutputProfiles(:,6) = interp1(Depths,ProfilesF(:,2),OutputDepths,'spline');
    
    MinOutput = MinimisationFunction(SearchVarFinal);
    
    %%% Minimisation subfunction
    function Value = MinimisationFunction(SearchVar,varargin)
        
        % Get profiles
        [~,Profiles,InterfaceIndicies] = GetProfiles(IntDomFcn(SearchVar),InitCondFcn(SearchVar));
        
        % Calculate residual
        Value = ResidualCond(Profiles,InterfaceIndicies);
    end
    
    
    %%% Full domain solve
    function [Grid,Profiles,InterfaceIndicies] = GetProfiles(Domains,SurfaceProfiles)
        
        % Integrate through mush
        [Mushy_Grid,Mushy_Profiles] = MushyDiffusionEquations(SurfaceProfiles,Domains{1});
        
        % Integrate temperature through liquid
        [Liquid_Grid,Liquid_Temp] = SimpleDiffusionEquation(1,InterfaceCondTemp(Mushy_Profiles),Domains{2});
        
        % Combine profiles
        Grid = [Mushy_Grid;Liquid_Grid];
        Profiles = [Mushy_Profiles;Liquid_Temp,ones(size(Liquid_Grid))];
        InterfaceIndicies = size(Mushy_Profiles,1) + [0,1];
    end
    
    
    %%% Mushy diffusion integration subfunction
    function [Grid,Profiles] = MushyDiffusionEquations(InitialConditions,Domain)
        
        % Set diffusion equation
        ForcingTerms = @(z,Var) [1;-z/2;Var(3)]*Var(2);
        MassTerms = @(z,Var) [1,0,0; 0,1,z*Stefan/2;0,0,(LqdusGrad+LqdusInf-Var(1))];
        
        % Perform integration
        w = warning('off','all');
        [Grid,Profiles] = ode15s(ForcingTerms,Domain,InitialConditions,odeset('Mass',MassTerms));
        warning(w);
    end 
    
    
    %%% Simple diffusion integration subfunction
    function [Z,Prof] = SimpleDiffusionEquation(Diffusivity,InitialConditions,Domain)
        
        % Set diffusion equation
        DiffEqn = @(z,Var) [Var(2);-(z*Var(2))/(2*Diffusivity)];
        
        % Perform integration
        [Z,Prof] = ode45(DiffEqn,Domain,InitialConditions);
    end
    
    
    %%% Parameter unpacking subfunction
    function [Stefan,LqdusGrad,LqdusInf,DomainDepth] = ParamUnpack(Param,Comp)
        Stefan    = Param.Stefan;
        LqdusInf  = Param.LqdusInf;
        LqdusGrad = Param.LqdusGrad;
        DomainDepth = Comp.DomainDepth;
    end
end

