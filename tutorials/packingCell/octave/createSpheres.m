%%createSpheres
% Federico Municchi, Nottingham (2019)
% Script to create bcc, fcc or sic sphere configurations for OpenFOAM
% (blockMesh and snappyHexMesh)
% Always 3D!
close all;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT DATA
run('input.m');
config = 'fcc'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START COMPUTATION

%Calculate sphere diameter based on configuration
%Each confuguration has a different number of (full) particles in the cell
%See Fig.1 in 10.1103/PhysRevE.94.053118

nparFull=0;
if strcmp(config ,'fcc')
  display('This is FCC');
  nparFull = 4.0;
elseif strcmp(config ,'bcc')
  display('This is BCC');
  nparFull = 3.0;
else
  display('This is SC');
  nparFull = 1.0;
endif

%Particle diameter
d = power(phi*(L^3)*6.0/(nparFull*pi),1.0/3.0);
filename = "particleData.H";
fid = fopen (filename, "w");
fprintf(fid,"dp %f;\n",d);
fclose(fid);
%Now create vectors with particle centres matrix

if strcmp(config ,'fcc')
  %There are 14 particles in total
  pCentre(:,1) = [-1 -1 -1];
  pCentre(:,2) = [-1 1 -1];
  pCentre(:,3) = [-1 1 1];
  pCentre(:,4) = [-1 -1 1];
  pCentre(:,5) = [-1 0 0];
  pCentre(:,6) = [0 -1 0];
  pCentre(:,7) = [0  1 0];
  pCentre(:,8) = [0 0 -1];
  pCentre(:,9) = [0 0  1];

  for i=1:5
    pCentre(1,9+i) = pCentre(1,i) + 2.0 ;
    pCentre(2,9+i) = pCentre(2,i);
    pCentre(3,9+i) = pCentre(3,i);
  endfor
elseif strcmp(config ,'bcc')
  %There are 9 particles in total
  pCentre(:,1) = [-1 -1 -1];
  pCentre(:,2) = [-1 1 -1];
  pCentre(:,3) = [-1 1 1];
  pCentre(:,4) = [-1 -1 1];
  pCentre(:,5) = [0 0 0];

  for i=1:4
    pCentre(1,5+i) = pCentre(1,i) + 2.0 ;
    pCentre(2,5+i) = pCentre(2,i);
    pCentre(3,5+i) = pCentre(3,i);
  endfor
else
  %Only one particle in sic
  pCentre(:,1) = [0 0 0];
endif

%Scale with the cell size
pCentre = pCentre * L/2.0;

%Finally write first cell
filename = "spheres.H";
fid = fopen (filename, "w");
sphereId = 0;
for i=1:length(pCentre(1,:))
  writeSphereSnappy(fid,pCentre(:,i),d/2.0,sphereId);
  sphereId=sphereId+1;
endfor
fclose (fid);
