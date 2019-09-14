%% reactiveStationary.m -- an executable m-file for solving a boundary-value problem
% Automatically created in CHEBGUI by user pmzfm.
% Created on July 16, 2019 at 15:50.

%% Problem description.
% Solving
%   1.23*u'' - 0.0999*u' -2.79952* u = 0,
% for x in [0, 10], subject to
%   u(0) = 1,
%   u'(10) = 0.
clear all;
close all;

%% Problem set-up.
% Define the domain.
dom = [0, 16];
params = [  0.2997    0.1107    0.1278    0.9955; ...
            0.7131    0.5833    0.9340    0.9921; ...
            0.030183999	0.029070882	0.042879496	0.9997321; ...
            5.9140789	0.084172359	0.16638671	0.9081770; ...
            0.72843623	0.61017891	0.97419645	0.99190715 ...
          ];
porosity = 1; %It should always be =1 since the Favre velocity is 1 already

FileList = ["fcc_Pe10_Da_1.04.csv"; ...
            "fcc_Pe100_Da_962.csv"; ...
            "fcc_Pe100_Da_1.04.csv"; ...
            "fcc_Pe10_Da_962.csv"; ...
            "fcc_Pe100_Da_100000.csv" ...
            ];

caseId =5;

fileName = FileList(caseId);
M = csvread(fileName);


% Assign the differential equation to a chebop on that domain.
N = chebop(@(x,u) (params(caseId,2))*diff(u,2)-(porosity*params(caseId,4))*diff(u)-(params(caseId,1))*u, dom);

% Set up the rhs of the differential equation so that N(u) = rhs.
rhs = 0;

% Assign boundary conditions to the chebop.
N.bc = @(x,u) [u(dom(1))-1; feval(diff(u),dom(2))];

% Construct a linear chebfun on the domain, 
x = chebfun(@(x) x, dom);
% and assign an initial guess to the chebop.
N.init = 1;

%% Setup preferences for solving the problem.
% Create a CHEBOPPREF object for passing preferences.
% (See 'help cheboppref' for more possible options.)
options = cheboppref();

% Print information to the command window while solving:
options.display = 'iter';

% Option for tolerance.
options.bvpTol = 5e-13;

% Option for damping.
options.damping = false;

% Specify the discretization to use. Possible options are:
%  'values' (default)
%  'coeffs'
%  A function handle (see 'help cheboppref' for details).
options.discretization = 'values';

% Option for determining how long each Newton step is shown.
options.plotting = 0.1;

%% Solve!
% Call solvebvp to solve the problem.
% (With the default options, this is equivalent to u = N\rhs.)
u = solvebvp(N, rhs, options);

%% Plot the solution.
figure(1)
%M = csvread(fileName)
x = linspace(0 ,15,16);
plot(u, 'LineWidth',2 , '-k')
hold on;
plot (x,M(10:25,1)/M(10,1),'.k','markerSize',20)
l =legend('Upscaled model', ...
       'Fully resolved simulation', ...
        'location','Best');
legend boxoff;
set(gca,'FontSize',26);
set(l, 'interpreter', 'latex');
box on;
%set(gca, 'XScale', 'log');

%set(gca, 'YScale', 'log');
xlabel('$x$','interpreter','latex');
ylabel('$<c>$','interpreter','latex');
xlim([0 12]);
dim1 = [.6 .4 .2 .2];
dim2 = [.6 .5 .2 .2];
dim3 = [.6 .6 .2 .2];
annotation('textbox',dim1,'String','$\epsilon = 0.7 $','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman','interpreter','latex')
annotation('textbox',dim2,'String','$\mathrm{Pe}=100$','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman','interpreter','latex')
annotation('textbox',dim3,'String','$\mathrm{Da}_{\mathrm{II}} = 100000$','FitBoxToText','on','EdgeColor',[1 1 1,],'FontSize',30, 'FontName', 'Times New Roman','interpreter','latex')

set(gca, 'FontName', 'Times New Roman');

x0=10;
y0=10;
width=1100;
height=800;
set(gcf,'position',[x0,y0,width,height])
saveas(gcf,'validation','epsc')
hold off;
