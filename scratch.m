%%
runOpt_weight;
runOpt_severity;
runOpt_resolution;

%% 3/4/2015 Combine four .mat files with 10 noise realizations each 
% It would be nice to 40 noise realizations in one run. However, it will run out of memory. 
clear;close all;
load('severity_multiNoise1.mat');
pVals_s1=pVals_s;
fVals_s1=fVals_s;
dscL_s1=dscL_s;
measuredImgs_s1=measuredImgs_s;
figure;plot(dscL_s1(1,:),'r*-');

load('severity_multiNoise2.mat');
pVals_s2=pVals_s;
fVals_s2=fVals_s;
dscL_s2=dscL_s;
measuredImgs_s2=measuredImgs_s;
hold on;plot(dscL_s2(1,:),'m*-');

load('severity_multiNoise3.mat');
pVals_s3=pVals_s;
fVals_s3=fVals_s;
dscL_s3=dscL_s;
measuredImgs_s3=measuredImgs_s;
hold on;plot(dscL_s3(1,:),'k*-');

load('severity_multiNoise4.mat');
pVals_s4=pVals_s;
fVals_s4=fVals_s;
dscL_s4=dscL_s;
measuredImgs_s4=measuredImgs_s;
hold on;plot(dscL_s4(1,:),'b*-');

% The plot confirms the difference between noise realizations

nNoise_s=40;
ns=6;
pVals_s=zeros(nNoise_s,length(tp_s),ns);
fVals_s=zeros(ns,nNoise_s);
dscL_s=zeros(ns,nNoise_s);
measuredImgs_s=zeros(dimX,dimY,ns,nNoise_s);

pVals_s(1:10,:,:)=pVals_s1;
pVals_s(11:20,:,:)=pVals_s2;
pVals_s(21:30,:,:)=pVals_s3;
pVals_s(31:40,:,:)=pVals_s4;

fVals_s(:,1:10)=fVals_s1;
fVals_s(:,11:20)=fVals_s2;
fVals_s(:,21:30)=fVals_s3;
fVals_s(:,31:40)=fVals_s4;

dscL_s(:,1:10)=dscL_s1;
dscL_s(:,11:20)=dscL_s2;
dscL_s(:,21:30)=dscL_s3;
dscL_s(:,31:40)=dscL_s4;

measuredImgs_s(:,:,:,1:10)=measuredImgs_s1;
measuredImgs_s(:,:,:,11:20)=measuredImgs_s2;
measuredImgs_s(:,:,:,21:30)=measuredImgs_s3;
measuredImgs_s(:,:,:,31:40)=measuredImgs_s4;

clear dscL_s1 dscL_s2 dscL_s3 dscL_s4
clear fVals_s1 fVals_s2 fVals_s3 fVals_s4
clear pVals_s1 pVals_s2 pVals_s3 pVals_s4
clear measuredImgs_s1 measuredImgs_s2 measuredImgs_s3 measuredImgs_s4
save('severity_multiNoise.mat');


