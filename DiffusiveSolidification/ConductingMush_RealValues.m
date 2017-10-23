%% Real Values Investigation
%
% (14/01/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Setup Investigation

% Settings
DeltaTempVec = [5,10,20];
StefanStart = 10.^-2;
StefanFactor = 10.^0.2;
DepthScale = 0.01; % Meters
Times = [1,10,60]; % Days
ThermDiff = 1.3e-7; % Meters ^2 / Second

% Computational Domain
Comp.DomainDepth = 15;
Comp.NGrid       = [59,49];

%% Get results

% Prepare containers
GrowthRate = NaN(1,length(DeltaTempVec)+1);
Grid = NaN(Comp.NGrid*[1;1]+1,length(DeltaTempVec)+1);
LqdFrac = NaN(Comp.NGrid*[1;1]+1,length(DeltaTempVec)+1);

% Loop over temperature ranges
for iD = 1:length(DeltaTempVec)+1
    
    % Set variables
    if iD <= length(DeltaTempVec)
        LqdusInf = 1-2/DeltaTempVec(iD);
        LqdusGrad = 3/DeltaTempVec(iD);
        StefanTarget = 83.4/DeltaTempVec(iD);
    else % Wettlaufer comparison
        LqdusInf = (-4 - -20)/(-2 - -20);
        LqdusGrad = 0.085*70/(-2 - -20);
        StefanTarget = 83.4/(-2 - -20);
    end
    
    % Loop over Stefan numbers
    Stefan = StefanStart;
    IG  = [DiffGrow_Simp_Cond(1,LqdusInf),1/sqrt(pi),LqdusGrad/(LqdusGrad+LqdusInf)];
    while 1
        
        % Set variables
        Stefan = min(StefanTarget,Stefan*StefanFactor);
        IG(3) = LqdusGrad/(LqdusGrad+LqdusInf);
        
        % Get results
        [GrowthRate(iD),Profiles,Grid(:,iD),IG] = DiffGrow_Mushy_Cond_NS(Stefan,LqdusGrad,LqdusInf,Comp,IG);
        LqdFrac(:,iD) = Profiles(:,5);
        
        % Exit loop
        if Stefan == StefanTarget
            disp('Loop Exit')
            break
        end
    end
    
    
end

% Tidy
clear iD StefanTarget Stefan IG Profiles


%% Process results

% Remove non-mush
Grid = Grid(1:Comp.NGrid(1)+1,:);
LqdFrac = LqdFrac(1:Comp.NGrid(1)+1,:);

% Find half-solid depth
HalfLiquidFraction = NaN(1,3);
for iD = 1:length(DeltaTempVec)
    HalfLiquidFraction(iD) = 1-interp1(LqdFrac(:,iD),Grid(:,iD)/GrowthRate(iD),0.5);
end

% Find physical depths
PhysDpths = sqrt(ThermDiff*24*60^2*Times')*GrowthRate/DepthScale;

