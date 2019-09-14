clear all;
close all;

%%============================================================================%%
%% INPUT

run('ranges.m')
config = 'fcc';
generateStudy = true;
%%============================================================================%%
%% CREATE CASES
root = pwd;

if generateStudy
    command = sprintf('rm -r %s_*',config);
    system(command);
    for n=1:length(phiRange)

        caseName = sprintf('%s_%.1f',config,phiRange(n));
        command = sprintf('cp -r ../../cases/specCellFoam/packingCell %s',caseName);
        system(command);
        cd(caseName)

        fid = fopen('./octave/input.m','w')
        fprintf(fid,'phi = %f;\n',phiRange(n));
        fprintf(fid,'L = 1.0;\n');
        fprintf(fid,'NCells =1;\n');
        fclose(fid);

        system('./createMesh.sh');
        system('./runFlow.sh');

        cd(root);

    endfor
endif

%%============================================================================%%
%% RUN PARAMETRIC STUDY

% Create matrix phi Pe Da lambda Dx Ux
M = zeros(length(phiRange)*length(PeRange)*length(DaRange),7);
% Create same matrix for Dirichlet BCs
D = zeros(length(phiRange)*length(PeRange),6);
% Create same matrix for Neumann BCs
N = zeros(length(phiRange)*length(PeRange),6);
id=1;
idD=1;
for n=1:length(phiRange)

    caseName = sprintf('%s_%.1f',config,phiRange(n));
    cd(caseName)
    for pn=1:length(PeRange)

        %First run Dirichlet
        fid = fopen('./constant/DA.H','w')
        fprintf(fid,'Da  -1;\n');
        fprintf(fid,'Pe  %f;\n',PeRange(pn));
        fprintf(fid,'BC  fixedValue;\n');
        fprintf(fid,'val 1;\n');
        fclose(fid);
        command = sprintf('rm log.specCellFoam specCellFoam_Diric_Pe%i.log',PeRange(pn));
        system(command);
        system('./runSpecFoam');

        X = csvread('results.csv');
        D(idD,:) = [phiRange(n) PeRange(pn)  X(1,1) X(1,2) X(1,3) X(1,4)];

        
        command = sprintf('mv log.specCellFoam specCellFoam_Diric_Pe%i.log',PeRange(pn));
        system(command);
        
        %Now run Neumann
        fid = fopen('./constant/DA.H','w')
        fprintf(fid,'Da  -1;\n');
        fprintf(fid,'Pe  %f;\n',PeRange(pn));
        fprintf(fid,'BC  zeroGradient;\n');
        fprintf(fid,'val 1;\n');
        fclose(fid);
        command = sprintf('rm log.specCellFoam specCellFoam_Neu_Pe%i.log',PeRange(pn));
        system(command);
        
        system('./runTransportUpscaling');

        command = sprintf('mv log.specCellFoam specCellFoam_Neu_Pe%i.log',PeRange(pn));
        system(command);

        X = csvread('results.csv');
        N(idD,:) = [phiRange(n) PeRange(pn)  X(1,1) X(1,2) X(1,3) X(1,4)];


        idD=idD+1;

        for dn=1:length(DaRange)

            fid = fopen('./constant/DA.H','w')
            fprintf(fid,'Da  %f;\n',-DaRange(dn));
            fprintf(fid,'Pe  %f;\n',PeRange(pn));
            fprintf(fid,'BC Robin;\n');
            fprintf(fid,'val 1;\n');
            fclose(fid);
            command = sprintf('rm log.specCellFoam specCellFoam_Da%i_Pe%i.log',DaRange(dn),PeRange(pn));
            system(command);
            system('./runSpecFoam');
            command = sprintf('mv log.specCellFoam specCellFoam_Da%i_Pe%i.log',DaRange(dn),PeRange(pn));
            system(command);

            X = csvread('results.csv');

            M(id,:) = [phiRange(n) PeRange(pn) -DaRange(dn) X(1,1) X(1,2) X(1,3) X(1,4)];
            id=id+1;

        endfor
    endfor


    cd(root);

endfor

%Save to file
save  -mat7-binary final_results.mat M;
save  -mat7-binary final_Dirichlet.mat D;
save  -mat7-binary final_Neumann.mat N;
