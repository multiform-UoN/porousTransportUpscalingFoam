clear all;
close all;
lineWidth = 2;
load('final_results.mat','-mat');
load('final_Dirichlet2.mat','-mat');
load('final_Neumann.mat','-mat');
run('ranges.m')
set(0,'defaultTextInterpreter','latex'); %trying to set the default
markerSize = 15;
stride = length(phiRange)*length(PeRange)*length(DaRange);

phiId = 0;
stride = phiId*length(PeRange)*length(DaRange)+1;
strideN = phiId*length(PeRange)+1;
DaN = length(DaRange);
hold on;
figure(1)
semilogx(DaRange,(M((stride):(stride+DaN-1),4)-N(strideN,3))/(D(strideN,3)-N(strideN,3)),'-k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+DaN):(stride+2*DaN-1),4)-N(strideN+1,3))/(D(strideN+1,3)-N(strideN+1,3)),'--k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+2*DaN):(stride+3*DaN-1),4)-N(strideN+2,3))/(D(strideN+2,3)-N(strideN+2,3)),'-.k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+3*DaN):(stride+4*DaN-1),4)-N(strideN+3,3))/(D(strideN+3,3)-N(strideN+3,3)),':k','linewidth',lineWidth);
%semilogx(DaRange,(M((stride+4*DaN):(stride+5*DaN-1),4)-N(strideN+4,3))/(D(strideN+4,3)-N(strideN+4,3)),'-+k','linewidth',lineWidth,'markersize',markerSize);
l =legend('$\mathrm{Pe} = 0.1$', ...
       '$\mathrm{Pe} = 1$', ...
       '$\mathrm{Pe} = 10$', ...
       '$\mathrm{Pe} = 100$', ...    
        'location','Best');
legend boxoff;
set(gca,'FontSize',26);
set(l, 'interpreter', 'latex');
box on;
set(gca, 'XScale', 'log');

%set(gca, 'YScale', 'log');
xlabel('$\mathrm{Da}_{\mathrm{II}}$','interpreter','latex');
ylabel('$\lambda^{\prime}$','interpreter','latex');
xlim([0.01 1e+05]);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String','\epsilon = 0.9','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman');
x0=10;
y0=10;
width=1100;
height=800;
set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'lambda','epsc')
hold off;



figure(2)
semilogx(DaRange,(M((stride):(stride+DaN-1),5)-N(strideN,4))/(D(strideN,4)-N(strideN,4)),'-k','linewidth',lineWidth);

hold on;
semilogx(DaRange,(M((stride+DaN):(stride+2*DaN-1),5)-N(strideN+1,4))/(D(strideN+1,4)-N(strideN+1,4)),'--k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+2*DaN):(stride+3*DaN-1),5)-N(strideN+2,4))/(D(strideN+2,4)-N(strideN+2,4)),'-.k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+3*DaN):(stride+4*DaN-1),5)-N(strideN+3,4))/(D(strideN+3,4)-N(strideN+3,4)),':k','linewidth',lineWidth);
%semilogx(DaRange,(M((stride+4*DaN):(stride+5*DaN-1),4)-N(strideN+4,3))/(D(strideN+4,3)-N(strideN+4,3)),'-+k','linewidth',lineWidth,'markersize',markerSize);
l =legend('$\mathrm{Pe} = 0.1$', ...
       '$\mathrm{Pe} = 1$', ...
       '$\mathrm{Pe} = 10$', ...
       '$\mathrm{Pe} = 100$', ...    
        'location','Best');
legend boxoff;
set(gca,'FontSize',26);
set(l, 'interpreter', 'latex');
box on;
set(gca, 'XScale', 'log');

%set(gca, 'YScale', 'log');
xlabel('$\mathrm{Da}_{\mathrm{II}}$','interpreter','latex');
ylabel('$\mathcal{D}^{\prime}_{xx}$','interpreter','latex');
xlim([0.01 1e+05]);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String','\epsilon = 0.9 ','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman');
x0=10;
y0=10;
width=1100;
height=800;
set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'Dxx','epsc')
hold off;

figure(3)
semilogx(DaRange,(M((stride):(stride+DaN-1),6)-N(strideN,5))/(D(strideN,5)-N(strideN,5)),'-k','linewidth',lineWidth);

hold on;
semilogx(DaRange,(M((stride+DaN):(stride+2*DaN-1),6)-N(strideN+1,5))/(D(strideN+1,5)-N(strideN+1,5)),'--k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+2*DaN):(stride+3*DaN-1),6)-N(strideN+2,5))/(D(strideN+2,5)-N(strideN+2,5)),'-.k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+3*DaN):(stride+4*DaN-1),6)-N(strideN+3,5))/(D(strideN+3,5)-N(strideN+3,5)),':k','linewidth',lineWidth);
%semilogx(DaRange,(M((stride+4*DaN):(stride+5*DaN-1),4)-N(strideN+4,3))/(D(strideN+4,3)-N(strideN+4,3)),'-+k','linewidth',lineWidth,'markersize',markerSize);
l =legend('$\mathrm{Pe} = 0.1$', ...
       '$\mathrm{Pe} = 1$', ...
       '$\mathrm{Pe} = 10$', ...
       '$\mathrm{Pe} = 100$', ...    
        'location','Best');
legend boxoff;
set(gca,'FontSize',26);
set(l, 'interpreter', 'latex');
box on;
set(gca, 'XScale', 'log');

%set(gca, 'YScale', 'log');
xlabel('$\mathrm{Da}_{\mathrm{II}}$','interpreter','latex');
ylabel('$\mathcal{D}^{\prime}_{yy}$','interpreter','latex');
xlim([0.01 1e+05]);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String','\epsilon = 0.9 ','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman');
x0=10;
y0=10;
width=1100;
height=800;
set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'Dyy','epsc')
hold off;

figure(4)
semilogx(DaRange,(M((stride):(stride+DaN-1),7)-N(strideN,6))/(D(strideN,6)-N(strideN,6)),'-k','linewidth',lineWidth);
hold on;
semilogx(DaRange,(M((stride+DaN):(stride+2*DaN-1),7)-N(strideN+1,6))/(D(strideN+1,6)-N(strideN+1,6)),'--k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+2*DaN):(stride+3*DaN-1),7)-N(strideN+2,6))/(D(strideN+2,6)-N(strideN+2,6)),'-.k','linewidth',lineWidth);
semilogx(DaRange,(M((stride+3*DaN):(stride+4*DaN-1),7)-N(strideN+3,6))/(D(strideN+3,6)-N(strideN+3,6)),':k','linewidth',lineWidth);
%semilogx(DaRange,(M((stride+4*DaN):(stride+5*DaN-1),4)-N(strideN+4,3))/(D(strideN+4,3)-N(strideN+4,3)),'-+k','linewidth',lineWidth,'markersize',markerSize);
l =legend('$\mathrm{Pe} = 0.1$', ...
       '$\mathrm{Pe} = 1$', ...
       '$\mathrm{Pe} = 10$', ...
       '$\mathrm{Pe} = 100$', ...    
        'location','Best');
legend boxoff;
set(gca,'FontSize',26);
set(l, 'interpreter', 'latex');
box on;
set(gca, 'XScale', 'log');

%set(gca, 'YScale', 'log');
xlabel('$\mathrm{Da}_{\mathrm{II}}$','interpreter','latex');
ylabel('$V^{\prime}_{x}$','interpreter','latex');
xlim([0.01 1e+05]);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String','\epsilon = 0.9 ','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman');
x0=10;
y0=10;
width=1100;
height=800;
set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'Vx','epsc')
hold off;