function [InterfacePosition,varargout] = DiffGrow_Simp_Biot_NS(Stefan,LqdusGrad,LqdusInf,Biot,varargin)
% This function solves for the interface position in the simplified model
% for a Biot number dependent boundary cooling but without salt diffusion
% by solving the implicit equation directly.
%
% Inputs: MushDiff  - Mushy diffusivity factor
%         LqdusGrad - Liquidus gradient
%         LqdusInf  - Initial liquidus temperature
%         Biot      - Biot number
% (Opt.)  Time      - Time
%
% (Opt.)  Comp  - Computation structure
%                    Required fields:
%                      DomainDepth - Computational domain depth
%                      NGrid       - Number of grid points in each domain
%                                    (Two values)
% 
% Outputs: InterfacePosition - Scaled position of the mush/liquid
%                                interface.
%          Depths            - Depths associated with the output profiles
%          OutputProfiles    - Enthalpy, temperature, bulk salinity, liquid
%                                salinity and liquid fraction profiles.
%
% (17/03/15)

    % Set mushy diffusivity
    MushDiff     = 1/(1+Stefan/LqdusGrad);
    MushDiffFact = sqrt(1/MushDiff);

    % Set time
    if nargin > 4
        Time = varargin{1};
    else
        Time = 1;
    end
    
    % Set 'Comp' input
    CalcProfs = false;
    if nargin > 5
        CalcProfs = true;
        Comp = varargin{2};
    elseif nargout > 1
        error('Must supply ''Comp'' argument to get profiles')
    end
    
    % Check for freezing
    if LqdusInf >= BiotCoolingProfile(0,Time,Biot,1)
        
        % Imperfect cooling
        if ~isinf(Biot)
            
            % Get self-similar Biot number
            ssBiot = Biot*sqrt(Time);
            
            % Set interface condition function
            erfx = @(y) erf(y).*exp(y.^2);
            IntCondFcn = @(Pos) LqdusInf*(erfcx(Pos/2)/erfcx(Pos/2+ssBiot)-1) ...        
                - (1-LqdusInf)*(erfx(MushDiffFact*Pos/2)/erfcx(MushDiffFact*Pos/2+ssBiot/MushDiffFact)+1);
            
            % Find interface
            InterfacePosition = fzero(IntCondFcn,0.01);
            
        % Perfect cooling
        else
            
            % Call routine
            InterfacePosition = DiffGrow_Simp_Cond(MushDiffFact,LqdusInf);
        end
        
        % Invert self-similarity transformation
        InterfacePosition = InterfacePosition*sqrt(Time);
        
        % Set profiles output
        if CalcProfs
            
            % Set domain depths
            Mushy_Depths  = ChebGrid(Comp.NGrid(1),InterfacePosition);
            Liquid_Depths = InterfacePosition + ChebGrid(Comp.NGrid(2),Comp.DomainDepth-InterfacePosition);
            Depths        = [Mushy_Depths;Liquid_Depths(2:end)];
            
            % Set profiles
            OutputProfiles = NaN(length(Depths),5);
            OutputProfiles(1:Comp.NGrid(1)+1,2)   = LqdusInf * BiotCoolingProfile(Mushy_Depths,Time,Biot,MushDiff) / BiotCoolingProfile(InterfacePosition,Time,Biot,MushDiff);
            OutputProfiles(Comp.NGrid(1)+2:end,2) = 1-(1-LqdusInf)*(1-BiotCoolingProfile(Liquid_Depths(2:end),Time,Biot,1))/(1-BiotCoolingProfile(InterfacePosition,Time,Biot,1));
            OutputProfiles(:,3) = ones(length(Depths),1);
            OutputProfiles(:,4) = max((LqdusInf-OutputProfiles(:,2))/LqdusGrad,0);
            OutputProfiles(:,5) = OutputProfiles(:,3)./(1+OutputProfiles(:,4));
            OutputProfiles(:,1) = OutputProfiles(:,2) - Stefan*(1-OutputProfiles(:,5));
        end
        
    % No freezing
    else
        
        % No interface
        InterfacePosition = NaN;
        
        % Set output profiles
        if CalcProfs
            
            % Set domain depths
            Mushy_Depths  = NaN(Comp.NGrid(1),1);
            Liquid_Depths = ChebGrid(Comp.NGrid(2),Comp.DomainDepth);
            Depths        = [Mushy_Depths;Liquid_Depths];
    
            % Pure thermal cooling
            OutputProfiles = NaN(length(Depths),5);
            OutputProfiles(1:Comp.NGrid(1),2)   = NaN(size(Mushy_Depths));
            OutputProfiles(Comp.NGrid(1)+1:end,2) = BiotCoolingProfile(Liquid_Depths,Time,Biot,1);
            OutputProfiles(:,3) = ones(length(Depths),1);
            OutputProfiles(:,4) = max((LqdusInf-OutputProfiles(:,2))/LqdusGrad,0);
            OutputProfiles(:,5) = OutputProfiles(:,3)./(1+OutputProfiles(:,4));
            OutputProfiles(:,1) = OutputProfiles(:,2) - Stefan*(1-OutputProfiles(:,5));
        end
    end
    
    % Set output
    if CalcProfs
        varargout{1} = Depths;
        varargout{2} = OutputProfiles;
    end
end

