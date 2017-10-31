%%
%
% This investigates the mushy layer growth in the Stefan regime.
% Note the calculation is repeated when studying mushy convection,
%
% It creates figure 3.17 of my pre-corrections thesis.
%
% (23/10/17)


clear
clc
addpath(genpath('../CoreFunctions'))


%% Parameters and settings

% Parameters
Stefan    = 300;
LqdusGrad = 10^-2;
LqdusInf  = 0.6;

% Computational settings
PDE.StartBiot = 10^-0.5;
PDE.MushyDynamics = true;
PDE.Depth = 1;
if 1 % High Res
    PDE.NGrid = 600;
    PDE.RelTol = 10^-7;
    PDE.AbsTol = 10^-8;
else % Low Res
    PDE.NGrid = 300;
    PDE.RelTol = 10^-4;
    PDE.AbsTol = 10^-6;
end

% Simple evolution settings
Comp.DomainDepth = 10;
Comp.NGrid       = [59,50];


%% Main Code

% Calculate freezing Biot
InterfaceTempFcn = @(Biot) LqdusInf-BiotCoolingProfile(0,1,Biot,1);
FreezingBiot = fzero(InterfaceTempFcn,1);

% Call sweep routine
[BiotNumbers,IntPos,Depths,Temps,LqdFracs] = DiffGrow_Mushy_Biot_NS(Stefan,LqdusInf,LqdusGrad,PDE);


% Tidy
clear InterfaceTempFcn

%% Simple model

% Loop over Biot numbers
for iB = 1:length(BiotNumbers)
    
    % Calculate simplified growth rate
    [GrowthRates(iB),~,tmp] = DiffGrow_Simp_Biot_NS(Stefan,LqdusGrad,LqdusInf,BiotNumbers(iB),1,Comp); %#ok<PFBNS>
    
	% First freezing point
	if iB == 1
        GrowthRates(iB) = 0;
    end
end


%% Output
FC = 0;
LW = 3;
FS = 18;
FW = 'normal';
PSFrag = true;

% Liquid Fractions
if 1
    
    % Options
    XTicks = 10.^(-1:4);%[0.1,0.3,1,3,10,30,100,300,1000,3000];
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    DepthGrid   = (0:0.001:1);
    LqdFracData = NaN(length(DepthGrid),length(BiotNumbers));
    for iB = 1:length(BiotNumbers)
        LqdFracData(:,iB) = interp1(Depths(:,iB),LqdFracs(:,iB),DepthGrid,'pchip');
    end
    
    % Create graph
    hold all
    contourf(BiotNumbers,DepthGrid,   LqdFracData,0:0.1:1,'edgecolor','k')
    plot([FreezingBiot;BiotNumbers(~isnan(IntPos))],[0;IntPos(~isnan(IntPos))],'r-' ,'linewidth',LW)
    contour(BiotNumbers,DepthGrid,   LqdFracData,0.5,'edgecolor',[1,0.5,0],'linewidth',LW)
    plot(FreezingBiot*[1,1],[0,1],'--','color',0.5*[1,1,1],'linewidth',LW)
    plot(1000*[1,1],[0,1],'--','color',0.5*[1,1,1],'linewidth',LW)
    plot(BiotNumbers,GrowthRates,'r--','linewidth',LW)
    plot(XTicks([1,1,end,end,1]),[0,0.08,0.08,0,0],'k')
    hold off
    if ~PSFrag
        xlabel('Biot number,  B_i','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',XTicks([1,end]),'xtick',XTicks)
        ylabel('Self-similar depth,  z / t^{1/2}','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','reverse','yscale','linear','ylim',[0,0.08],'ytick',0:0.02:0.08)
        caxis([0,1])
        colormap(cbrewer('seq','Blues',20))
        hC = colorbar('ydir','reverse','ylim',[0,1],'ytick',0:0.2:1);
        ylabel(hC,'Liquid fraction,  \chi','fontsize',FS,'fontweight',FW)
        line('parent',hC,'ydata',0.5*[1,1],'xdata',[-5,5],'color',[1,0.5,0],'linewidth',LW)
        line('parent',hC,'ydata',0.995*[1,1],'xdata',[-5,5],'color','r','linewidth',LW)
        set(gca,'fontsize',FS,'fontweight',FW)
        box on
    else
        xlabel('XLAB','fontsize',FS,'fontweight',FW)
        set(gca,'xscale','log','xlim',XTicks([1,end]),'xtick',XTicks,'xticklabel',{'X1','X2','X3','X4','X5','X6'})
        ylabel('YLAB','fontsize',FS,'fontweight',FW)
        set(gca,'ydir','reverse','yscale','linear','ylim',[0,0.08],'yticklabel',{'Y1','Y2','Y3','Y4','Y5'})
        caxis([0,1])
        colormap(cbrewer('seq','Blues',20))
        hC = colorbar('ydir','reverse','ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'C1','C2','C3','C4','C5','C6'});
        ylabel(hC,'CLAB','fontsize',FS,'fontweight',FW)
        for iC = 0:0.1:1
            line('parent',hC,'ydata',iC*[1,1],'xdata',[-5,5],'color','k','linewidth',1)
        end
        line('parent',hC,'ydata',0.5*[1,1],'xdata',[-5,5],'color',[1,0.5,0],'linewidth',LW)
        line('parent',hC,'ydata',0.995*[1,1],'xdata',[-5,5],'color','r','linewidth',LW)
        set(gca,'fontsize',FS,'fontweight',FW)
        box on
    end
    
    % Tidy
    clear XTicks DepthGrid LqdFracData iB hC
    
elseif 1
    
    % Options
    XTicks = 10.^(-2:3);
    YTicks = 0:0.2:1.2;
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    DepthGrid   = (0:0.01:1.5);
    LqdFracData = NaN(length(DepthGrid),length(BiotNumbers));
    for iB = 1:length(BiotNumbers)
        LqdFracData(:,iB) = interp1(Depths(:,iB),LqdFracs(:,iB),DepthGrid,'pchip');
    end
    
    % Create graph
    hold all
    contourf(BiotNumbers,DepthGrid,   LqdFracData,0:0.1:1,'edgecolor','k')
    plot([FreezingBiot;BiotNumbers(~isnan(IntPos))],[0;IntPos(~isnan(IntPos))],'r-' ,'linewidth',LW)
    contour(BiotNumbers,DepthGrid,   LqdFracData,0.5,'edgecolor',[1,0.5,0],'linewidth',LW)
    plot(XTicks([1,1,end,end,1]),YTicks([1,end,end,1,1]),'k')
    hold off
    xlabel('Biot number,  B_i','fontsize',FS,'fontweight',FW)
    set(gca,'xscale','log','xlim',XTicks([1,end]),'xtick',XTicks)
    ylabel('Self-similar depth,  z / t^{1/2}','fontsize',FS,'fontweight',FW)
    set(gca,'ydir','reverse','yscale','linear','ylim',YTicks([1,end]),'ytick',YTicks)
    caxis([0,1])
    colormap(cbrewer('seq','Blues',10))
    hC = colorbar('ydir','reverse','ylim',[0,1],'ytick',0:0.2:1);
    ylabel(hC,'Liquid fraction,  \chi','fontsize',FS,'fontweight',FW)
    line('parent',hC,'ydata',0.5*[1,1],'xdata',[-5,5],'color',[1,0.5,0],'linewidth',LW)
    line('parent',hC,'ydata',0.995*[1,1],'xdata',[-5,5],'color','r','linewidth',LW)
    set(gca,'fontsize',FS,'fontweight',FW)
    box on
    
    Ax1 = gca;
    Ax2 = axes('position',get(Ax1,'position'),'color','none');
    set(Ax2,'xscale','log','xlim',XTicks([1,end]),'xtick',XTicks,'xticklabel',[])
    set(Ax2,'ydir','reverse','ylim',YTicks([1,end]),'ytick',YTicks,'yticklabel',[])
    box on
    
    % Tidy
    clear XTicks DepthGrid LqdFracData iB hC
end


% Heat Capacities
if 0

    % Options
    XTicks = [0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000];
    CTicks = [1,3,10,30,100,300,1000,3000,10000,30000];
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    DepthGrid   = (0:0.001:1);
    HeatCapData = NaN(length(DepthGrid),length(BiotNumbers));
    for iB = 1:length(BiotNumbers)
        HeatCapData(:,iB) = interp1(Depths(:,iB),LqdFracs(:,iB),DepthGrid,'pchip');
    end
    HeatCapData = 1+ Stefan/LqdusGrad*HeatCapData.^2;
    HeatCapData1 = HeatCapData;
    HeatCapData2 = HeatCapData;
    HeatCapData1(HeatCapData == 30001) = 1;
    HeatCapData2(HeatCapData == 30001) = NaN;
    
    % Create graph
    hold all
    contourf(BiotNumbers,DepthGrid,log10(HeatCapData1),20,'edgecolor','none')
    plot([FreezingBiot;BiotNumbers(~isnan(IntPos))],[0;IntPos(~isnan(IntPos))],'r-' ,'linewidth',LW)
    contour(BiotNumbers,DepthGrid,log10(HeatCapData2),[2.5,3.5],'edgecolor','k','linewidth',LW)
    hold off
    xlabel('Biot number,  B_i','fontsize',FS)
    set(gca,'xscale','log','xlim',XTicks([1,end]),'xtick',XTicks)
    ylabel('Depth,  \eta','fontsize',FS)
    set(gca,'ydir','reverse','yscale','linear','ylim',[0,0.08],'ytick',0:0.01:0.08)
    caxis([0,4.5])
    colormap(cbrewer('seq','YlGnBu',20))
    hC = colorbar('ydir','reverse','ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',CTicks);
    ylabel(hC,'Effective heat capacity,  C_p^{mush}','fontsize',FS)
    line('parent',hC,'ydata',2.5 *[1,1],'xdata',[-5,5],'color','k','linewidth',LW)
    line('parent',hC,'ydata',3.5 *[1,1],'xdata',[-5,5],'color','k','linewidth',LW)
    line('parent',hC,'ydata',4.46*[1,1],'xdata',[-5,5],'color','r','linewidth',LW)
    set(gca,'fontsize',FS)
    box on
    
    % Tidy
    clear XTicks CTicks DepthGrid HeatCapData HeatCapData1 HeatCapData2 iB hC
end

% Tidy
clear FC LW FS FW