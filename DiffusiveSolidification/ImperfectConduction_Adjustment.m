%% Boundary Effects Script
% This script investigates the boundary effects for the simplified model.
%
% It produces figures 3.14, 3.15 and 3.16 (a,b&c) of my pre-corrections
% thesis.
%
% (23/10/2017)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Settings

% Parameter settings
InvLqdusGrad     = 10;
InvLqdusInfVec   = [0.03:0.01:0.1,0.12:0.02:0.90,0.91:0.01:0.99];
InvMushyDiffVec  = 10.^(0:0.1:3);
InvBiSamples     = 100;
InvBiMax         = 2*10^4;

% Simple evolution settings
SELqdusInf         = 0.8;
SEStefan           = 100;
SEBiVec            = 10.^(-2:0.1:4);
SEComp.DomainDepth = 10;
SEComp.NGrid       = [39,30];


%% Run investigation

% Prepare result variables
InvFreezingBiot = NaN(size(InvLqdusInfVec'));
InvBiotArray    = NaN(InvBiSamples+1,                        length(InvLqdusInfVec));
InvGrowthRates  = NaN(InvBiSamples+1,length(InvMushyDiffVec),length(InvLqdusInfVec));
InvSurfTemp     = NaN(InvBiSamples+1,length(InvMushyDiffVec),length(InvLqdusInfVec));
InvExponent     = NaN(               length(InvMushyDiffVec),length(InvLqdusInfVec));

% Loop over initial liquidus temperatures
tic
for iLi = 1:length(InvLqdusInfVec);
    
    % Calculate freezing Biot number
    FreezingBiotFcn   = @(Biot) InvLqdusInfVec(iLi) - BiotCoolingProfile(0,1,Biot,1);
    InvFreezingBiot(iLi) = fzero(FreezingBiotFcn,1);
    
    % Set investigation Biot numbers
    InvBiotArray(:,iLi) = [10.^(linspace(log10(InvFreezingBiot(iLi)),log10(InvBiMax),InvBiSamples)),inf]';
    
    % Loop over mushy diffusivities
    for iM = 1:length(InvMushyDiffVec)
        
        % Set parallel variables
        PAR_GrowthRates = NaN(InvBiSamples+1,1);
        PAR_SurfTemp    = NaN(InvBiSamples+1,1);
        
        % Loop over Biot numbers
        for iB = 1:InvBiSamples+1
            
            % Calculate simplified growth rate
            [PAR_GrowthRates(iB),~,tmp] = DiffGrow_Simp_Biot_NS(InvLqdusGrad*(InvMushyDiffVec(iM)^2-1),InvLqdusGrad,InvLqdusInfVec(iLi),InvBiotArray(iB,iLi),1,SEComp); %#ok<PFBNS>
            PAR_SurfTemp(iB) = tmp(1,2);
            
            % First freezing point
            if iB == 1
                PAR_GrowthRates(iB) = 0;
            end
        end
        
        % Perform fit
        ScaledGrowthRates = PAR_GrowthRates/PAR_GrowthRates(end);
        DataInd = (ScaledGrowthRates <= 0.6);
        DataInd = DataInd & [false; true(length(DataInd)-1,1)];
        SelectedBiot = InvBiotArray(DataInd,iLi)/InvFreezingBiot(iLi);
        SelectedData = ScaledGrowthRates(DataInd);
        FitFcn = @(p) 1-SelectedBiot.^-p;
        ErrFcn = @(p) sqrt(sum(((FitFcn(p)-SelectedData)./SelectedData).^2)/length(SelectedData));
        InvExponent(iM,iLi) = fminsearch(ErrFcn,-1);
        tmp_ExpErr(iM,iLi) = ErrFcn(InvExponent(iM,iLi));
        %tmp_Err2(iM,iLi) = sqrt( sum( ((FitFcn(InvExponent(iM,iLi))-SelectedData)./SelectedData).^2 )/length(SelectedData) );
        
        % Parallel output
        InvGrowthRates(:,iM,iLi) = PAR_GrowthRates;
        InvSurfTemp(   :,iM,iLi) = PAR_SurfTemp;
    end
    toc
end

% Tidy
clear iLi FreezingBiotFcn iM PAR_GrowthRates iB ScaledGrowthRates DataInd SelectedBiot SelectedData FitFcn ErrFcn


%% Calculate Simple Model Evolution

% Calculate freezing Biot number
FreezingBiotFcn = @(Biot) SELqdusInf - BiotCoolingProfile(0,1,Biot,1);
SEFreezingBiot  = fzero(FreezingBiotFcn,1);

% Prepare results container
SEIntPos  = NaN(1,length(SEBiVec));
SEDepths  = NaN(SEComp.NGrid*[1;1]+1,length(SEBiVec));
SELqdFrac = NaN(SEComp.NGrid*[1;1]+1,length(SEBiVec));

% Loop over Biot numbers
disp('Performing evolution in simple model')
tic
parfor iB = 1:length(SEBiVec)
    
    [SEIntPos(iB),SEDepths(:,iB),tmp] = DiffGrow_Simp_Biot_NS(SEStefan,InvLqdusGrad,SELqdusInf,SEBiVec(iB),1,SEComp);
    SELqdFrac(:,iB) = tmp(:,5);
end
SEDepths(isnan(SEDepths)) = -0.1;
toc

% Tidy
clear FreezingBiotFcn iB tmp


%% Output
FC = 0;
LW = 3;
FS = 18;
FW = 'normal';
PSFrag = true;


% Freezing Biot numbers (central range)
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    hold all
    plot(BiotCoolingProfile(0,1,10.^(-3:0.1:3),1),10.^(-3:0.1:3),'linewidth',LW)
    plot(NaN,NaN)
    plot(InvLqdusInfVec,5.8*exp(-4*InvLqdusInfVec),'--','linewidth',LW) % Trend by eye
    hold off
    if ~PSFrag
        xlabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS,'fontweight',FW)
        set(gca,'xtick',0:0.1:1)
        ylabel('Freezing Biot number,  B_f','fontsize',FS,'fontweight',FW)
        set(gca,'yscale','log','ylim',10.^[-3,3],'ytick',10.^(-3:3))
        set(gca,'fontsize',FS,'fontweight',FW)
        legend('Analytical result','Approximation')
    else
        xlabel('XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xtick',0:0.2:1,'xticklabel',{'X1','X2','X3','X4','X5','X6'})
        ylabel('YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'yscale','log','ylim',10.^[-3,3],'ytick',10.^(-3:3),'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6','Y7'})
        set(gca,'fontsize',FS,'fontweight',FW)
        legend('LEG111111111','LEG2')
    end
    axis square
    box on
end


% Freezing Biot numbers (low LqdusInf, not neat)
if 0
   
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    hold all
    plot(BiotCoolingProfile(0,1,10.^(-3:0.1:3),1),10.^(-3:0.1:3),'linewidth',LW)
    plot(NaN,NaN)
    plot(InvLqdusInfVec,1./(InvLqdusInfVec*sqrt(pi)),'--','linewidth',LW) % Trend by eye
    hold off
    xlabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS,'fontweight',FW)
    set(gca,'xscale','log','xtick',0:0.1:1)
    ylabel('Freezing Biot number,  B_f','fontsize',FS,'fontweight',FW)
    set(gca,'yscale','log','ylim',10.^[-3,3],'ytick',10.^(-3:3))
    set(gca,'fontsize',FS,'fontweight',FW)
    legend('Analytical result','Approximation')
    box on
end


% Freezing Biot numbers (high LqdusInf, not neat)
if 0
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    hold all
    plot(1-BiotCoolingProfile(0,1,10.^(-3:0.1:3),1),10.^(-3:0.1:3),'linewidth',LW)
    plot(NaN,NaN)
    plot(1-InvLqdusInfVec,sqrt(pi/4)*(1-InvLqdusInfVec),'--','linewidth',LW) % Trend by eye
    hold off
    xlabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS,'fontweight',FW)
    set(gca,'xdir','reverse','xscale','log','xtick',0:0.1:1)
    ylabel('Freezing Biot number,  B_f','fontsize',FS,'fontweight',FW)
    set(gca,'yscale','log','ylim',10.^[-3,3],'ytick',10.^(-3:3))
    set(gca,'fontsize',FS,'fontweight',FW)
    legend('Analytical result','Approximation')
    box on
end


% Evolution graph
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    hold all
    contourf(ones(SEComp.NGrid*[1;1]+1,1)*SEBiVec(1:end),SEDepths(:,1:end),SELqdFrac(:,1:end),0.9:0.005:1)
    plot(SEFreezingBiot*[1,1],[0,1],'--','color',0.5*[1,1,1],'linewidth',LW)
    plot(100*[1,1],[0,1],'--','color',0.5*[1,1,1],'linewidth',LW)
    plot([SEFreezingBiot,SEBiVec(~isnan(SEIntPos))],[0,SEIntPos(~isnan(SEIntPos))],'r','linewidth',LW)
    hold off
    if ~PSFrag
        xlabel('Biot number,  B_i','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xtick',10.^(-2:4))
        ylabel('Depth,  \eta','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','reverse','ylim',[0,1],'ytick',0:0.1:1)
        
        colormap(cbrewer('seq','Blues',40))
        caxis([0.9,1])
        hC = colorbar('ydir','reverse','ylim',[0.9,1],'ytick',0.9:0.01:1,'fontsize',FS,'fontweight',FW);
        ylabel(hC,'Liquid fraction,  \chi','fontsize',FS,'fontweight',FW)
        line('parent',hC,'xdata',[-5,5],'ydata',0.9993*[1,1],'color','r','linewidth',LW)
    else
        xlabel('XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xtick',10.^(-2:4),'xticklabel',{'X1','X2','X3','X4','X5','X6','X7'})
        ylabel('YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','reverse','ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6'})
        
        colormap(cbrewer('seq','Blues',40))
        caxis([0.9,1])
        hC = colorbar('ydir','reverse','ylim',[0.9,1],'ytick',0.9:0.02:1,'yticklabel',{'C1','C2','C3','C4','C5','C6'},'fontsize',FS,'fontweight',FW);
        ylabel(hC,'CLAB','fontsize',FS,'fontweight',FW)
        for iC = 0.9:0.005:1
            line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1],'color','k','linewidth',1)
        end
        line('parent',hC,'xdata',[-5,5],'ydata',0.9993*[1,1],'color','r','linewidth',LW)
    end
    set(gca,'fontsize',FS,'fontweight',FW)
    box on
    
    % Tidy
    clear hC
end


% Adjustment region against mushy diffusivity factor
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    Samp = 38;
    GrowthRateBiot = InvBiotArray(:,Samp)./InvFreezingBiot(Samp);
    GrowthRateData = InvGrowthRates(:,:,Samp)./(ones(InvBiSamples+1,1)*InvGrowthRates(end,:,Samp));
    
    % Create graph
    contourf(GrowthRateBiot(1:end-1),InvMushyDiffVec,GrowthRateData(1:end-1,:)',0:0.05:1)
    if ~PSFrag
        xlabel('Scaled Biot number,  B_i / B_f','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',10.^[0,3],'xtick',10.^(-2:4))
        ylabel('Mushy diffusivity factor,  m','fontsize',FS,'fontweight',FW)
        set(gca,'yscale','log','ytick',10.^(0:3))
        set(gca,'fontsize',FS,'fontweight',FW)
        
        caxis([0,1])
        colormap(cbrewer('seq','YlGnBu',40))
        ylabel(...
            colorbar('ylim',[0,1],'ytick',0:0.1:1,'fontsize',FS,'fontweight',FW), ...
            'Scaled growth rate,  \lambda(B_i) / \lambda(\infty)','fontsize',FS,'fontweight',FW)
    else
        xlabel('1XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',10.^[0,3],'xtick',10.^(0:3),'xticklabel',{'1X1','1X2','1X3','1X4'})
        ylabel('1YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'yscale','log','ytick',10.^(0:3),'yticklabel',{'1Y1','1Y2','1Y3','1Y4','1Y5','1Y6'})
        set(gca,'fontsize',FS,'fontweight',FW)
        
        caxis([0,1])
        colormap(cbrewer('seq','YlGnBu',40))
        hC = colorbar('ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'1C1','1C2','1C3','1C4','1C5','1C6'},'fontsize',FS,'fontweight',FW);
        ylabel(hC,'1CLAB','fontsize',FS,'fontweight',FW)
        for iC = 0:0.05:1
            line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1])
        end
    end
    
    % Tidy
    clear Samp GrowthRateBiot GrowthRateData
end


% Adjustment region against initial liquidus temperature
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    Samp = 9;
    GrowthRateBiot = InvBiotArray./(ones(InvBiSamples+1,1)*InvFreezingBiot');
    GrowthRateData = reshape(InvGrowthRates(:,Samp,:),[InvBiSamples+1,length(InvLqdusInfVec)])./(ones(InvBiSamples+1,1)*reshape(InvGrowthRates(end,Samp,:),[1,length(InvLqdusInfVec)]));
    
    % Create graph
    contourf(GrowthRateBiot(1:end-1,:)',InvLqdusInfVec'*ones(1,InvBiSamples),GrowthRateData(1:end-1,:)',0:0.05:1)
    patch(10.^[0,0,4,4],[0,InvLqdusInfVec([1,1]),0],0.5*[1,1,1])
    patch(10.^[0,0,4,4],[1,InvLqdusInfVec([end,end]),1],0.5*[1,1,1])
    if ~PSFrag
        xlabel('Scaled Biot number,  B_i / B_f','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',10.^[0,3],'xtick',10.^(0:3))
        ylabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS,'fontweight',FW)
        set(gca,'ylim',[0,1])
        set(gca,'fontsize',FS,'fontweight',FW)
    
        caxis([0,1])
        colormap(cbrewer('seq','YlGnBu',40))
        ylabel(...
            colorbar('ylim',[0,1],'ytick',0:0.1:1,'fontsize',FS,'fontweight',FW), ...
            'Scaled growth rate,  \lambda(B_i) / \lambda(\infty)','fontsize',FS,'fontweight',FW)
    else
        xlabel('1XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',10.^[0,3],'xtick',10.^(0:3),'xticklabel',{'1X1','1X2','1X3','1X4'})
        ylabel('2YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'2Y1','2Y2','2Y3','2Y4','2Y5','2Y6'})
        set(gca,'fontsize',FS,'fontweight',FW)
    
        caxis([0,1])
        colormap(cbrewer('seq','YlGnBu',40))
        hC = colorbar('ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'1C1','1C2','1C3','1C4','1C5','1C6'},'fontsize',FS,'fontweight',FW);
        ylabel(hC,'1CLAB','fontsize',FS,'fontweight',FW)
        for iC = 0:0.05:1
            line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1])
        end
    end
    
    % Tidy
    clear Samp GrowthRateBiot GrowthRateData
end


% Power law exponent
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    contourf(InvMushyDiffVec,InvLqdusInfVec,InvExponent',0.2:0.05:1)
    patch(InvMushyDiffVec([1,1,end,end]),[0,InvLqdusInfVec([1,1]),0],0.5*[1,1,1])
    patch(InvMushyDiffVec([1,1,end,end]),[1,InvLqdusInfVec([end,end]),1],0.5*[1,1,1])
    if ~PSFrag
        xlabel('Mushy Diffusivity Factor,  m','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log')
        ylabel('Initial Liquidus Temperature,  \theta_{L\infty}','fontsize',FS,'fontweight',FW)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1)
            set(gca,'fontsize',FS,'fontweight',FW)
    
        caxis([0.2,1])
        colormap(cbrewer('seq','YlGnBu',32))
        ylabel(...
            colorbar('ylim',[0.2,1],'ytick',0.2:0.1:1,'fontsize',FS,'fontweight',FW), ...
            'Power law exponent,  p','fontsize',FS,'fontweight',FW)
    else
        xlabel('3XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xtick',10.^(0:3),'xticklabel',{'3X1','3X2','3X3','3X4'})
        ylabel('3YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'3Y1','3Y2','3Y3','3Y4','3Y5','3Y6'})
            set(gca,'fontsize',FS,'fontweight',FW)
    
        caxis([0.2,1])
        colormap(cbrewer('seq','YlGnBu',32))
        hC = colorbar('ylim',[0.2,1],'ytick',0.2:0.2:1,'yticklabel',{'3C1','3C2','3C3','3C4','3C5'},'fontsize',FS,'fontweight',FW);
        ylabel(hC,'3CLAB','fontsize',FS,'fontweight',FW)
        for iC = 0:0.05:1
            line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1])
        end
    end
end


% Not pretty - Adjustment law chech
if 0
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Get data
    alc_MushyDiffIndex = 31;
    alc_LqdusInfIndex = 1:4:57;
    
    alc_Biot = InvBiotArray(:,alc_LqdusInfIndex);
    alc_GrowthRates = reshape(InvGrowthRates(:,alc_MushyDiffIndex,alc_LqdusInfIndex),101,length(alc_LqdusInfIndex));
    
    alc_ModGR = alc_Biot;
    alc_ModGR = (ones(101,1)*InvFreezingBiot(alc_LqdusInfIndex)')./alc_ModGR;
    alc_ModGR = 1-alc_ModGR.^(ones(101,1)*InvExponent(alc_MushyDiffIndex,alc_LqdusInfIndex));
    alc_ModGR = (ones(101,1)*alc_GrowthRates(end,:)).*alc_ModGR;
    
    % Create graph
    hold all
    plot(alc_Biot,alc_GrowthRates,'b','linewidth',3)
    plot(alc_Biot,alc_ModGR,'r--','linewidth',3)
    hold off
    set(gca,'xscale','log')
end


 % Tidy
 clear FC LW FS FW
 
 