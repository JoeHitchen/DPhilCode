function [BiotNumbers,IntPos,Depths,Temps,LqdFracs] = DiffGrow_Mushy_Biot_NS(Stefan,LqdusInf,LqdusGrad,PDE)
% This function performs multiple PDE integration runs and uses the results
% to give the behaviour of the self-similar solution as the Biot number is
% varied.
%
% Inputs:  Stefan - Stefan number
%          LqdusGrad - Liquidus gradient
%          LqdusInf  - Initial liquidus temperature
%          Biot      - Biot number
%          PDE       - PDE structure
%                      Required fields:
%                        zGrid         - Vertical grid (Not implimented)
%                        tGrid         - Time grid (Not implimented)
%                        MushyDynamics - Mushy (true) or simplified (false)
%                        RelTol        - Relative tolerance of iteration
%                        AbsTol        - Absolute tolerance of iteration
%
% Outputs: BiotNumbers - Self-similar Biot numbers of the run
%          IntPos      - Self-similar interface positions
%          Depths      - Self-similar depths
%          Temps       - Temperature profiles
%          LqdFracs    - Liquid fractions
%
% (22/03/16)


    % Settings (not yet pulled out)
    PDE.zGrid = PDE.Depth*(0:PDE.NGrid).^2/PDE.NGrid^2;
    %PDE.tGrid = [0,10.^(-1.2:0.05:-1)];
    PDE.tGrid = [0,10.^(-2:0.05:-1)];

    % Set start Biot number
    if ~isfield(PDE,'StartBiot')
        InterfaceTempFcn = @(Biot) LqdusInf-BiotCoolingProfile(0,1,Biot,1);
        StartBiot = fzero(InterfaceTempFcn,1);
    else
        StartBiot = PDE.StartBiot;
    end
    
    % Set Biot numbers for runs
    CoolFacStep = log10(sqrt(PDE.tGrid(end)/PDE.tGrid(2)));
    CoolFacVec = 10.^((4-log10(PDE.tGrid(end))/2):-CoolFacStep:log10(StartBiot));
    CoolFacVec = fliplr(CoolFacVec);
    
    % Initialise results containers
    TER_IntPos     = NaN(length(PDE.tGrid),length(CoolFacVec));
    TER_Temps      = NaN(length(PDE.zGrid),length(PDE.tGrid),length(CoolFacVec));
    
    % Loop over cooling factors
    hTic = tic;
    parfor iC = 1:length(CoolFacVec)
        
        % Call function
        [TER_IntPos(:,iC),TER_Temps(:,:,iC)] = DiffGrow_PDE_Routine(Stefan,LqdusGrad,LqdusInf,CoolFacVec(iC),PDE);
        
        % Show time elapsed
        toc(hTic)
    end
    
    % Set Biot numbers and interface positions
    BiotNumbers = sqrt(PDE.tGrid)'*CoolFacVec;
    IntPos = TER_IntPos./(sqrt(PDE.tGrid)'*ones(1,length(CoolFacVec)));
    
    % Combine interface positions
    BiotNumbers = [BiotNumbers(2,1);reshape(BiotNumbers(3:end,:),numel(BiotNumbers(3:end,:)),1)];
    IntPos      = [     IntPos(2,1);reshape(     IntPos(3:end,:),numel(     IntPos(3:end,:)),1)];
    
    % Combine temperature profiles
    Depths = NaN(length(PDE.zGrid),length(BiotNumbers));
    Temps  = NaN(length(PDE.zGrid),length(BiotNumbers));
    for iB = 1:length(BiotNumbers)
        
        % Set indexes
        iC = floor((iB-2)/(length(PDE.tGrid)-2))+1; iC(iC<1) = 1;
        it = iB-(length(PDE.tGrid)-2)*(iC-1)+1;
        
        % Set depths and temperatures
        Depths(:,iB) = PDE.zGrid./sqrt(PDE.tGrid(it));
        Temps( :,iB) = TER_Temps(:,it,iC);
    end
    
    % Set liquid fractions
    LqdFracs = LqdusGrad./(LqdusGrad+LqdusInf-Temps);
    LqdFracs(LqdFracs>1) = 1;
    LqdFracs(LqdFracs<0) = 1;
end