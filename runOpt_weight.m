%% This study is to evaluate the weight of constraint 
% True Parameter Settings
% Activity Ratio: Blood Pool : Myocardium : Background = 15 :95: 5
% Radius=20 mm, Thickness=8mm 
% For 256*256 matix size with 0.4mm*0.4mm pixel size, 
% Radius=50p, Thickness=20p
% PET resolution FWHM=6mm, Sigma=6.0mm/(2sqrt(2ln2))=6.0/2.3548=2.5480mm=6.37p

%% Defect=50,  50% severity, 40 degree extent
tp_w=[128.5 128.5 15 5 ...
    50 20 50 20 50 20 50 20 50 20 50 20 50 20 50 20 ...
    95 95 95 95 95 95 95 50 50 95 95 95 95 95 95 95 95 95];

global gaussFilter;
gaussFilter= fspecial('gaussian', [29 29], 6.37);
global imgMd;
global nRad; nRad=8;
global rAng; rAng=2*pi/nRad;
global hrAng; hrAng=pi/nRad;
global nSeg;nSeg=floor((length(tp_w)-4-2*nRad));
global sAng;sAng=2*pi/nSeg;
global dimX;dimX=256;
global dimY;dimY=256;

img_w=createActImg2D(tp_w);
figure;imshow(img_w,[]);title('Truth(Weight)');

%%% Poisson noise %%%
nImg_w= double(1e+15*imnoise(img_w*1e-15,'poisson')); 
trueGaussFilter=fspecial('gaussian', [29 29], 6.37);
imgMd=imfilter(nImg_w,trueGaussFilter,'same');
figure;imshow(imgMd,[]);title('Measured(Weight)');

initP_w=[127 126 20 8 ...
    45 27 46 26 44 28 45 30 47 26 46 25 43 29 44 30 ...
    70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70 70];

%imgInit_w=createActImg2D(initP_w);figure;imshow(imgInit_w,[]);
options_w = optimoptions(@fmincon,...
    'Display','iter',...
    'Algorithm','interior-point',...
    'FinDiffType','central',...
    'FinDiffRelStep',0.001,...
    'MaxFunEvals',10000 ...
    );

global weight;
weights=[0 0.1 0.5 1 5 10];
nw=length(weights);
pVals_w=zeros(nw,length(tp_w));
fVals_w=zeros(nw,1);
dscL_w=zeros(nw,1);
dmL_w=zeros(dimX,dimY,nw);
for k=1:nw
    weight=weights(k);
    [pVals_w(k,:),fVals_w(k)] = fmincon(@objConFunc,initP_w,[],[],[],[],[],[],@noconstraint,options_w); 
    [dscL_w(k), dmL_w(:,:,k)] = calcDSC(pVals_w(k,:), tp_w);
end;

formats_w=['ro-';  'bo-'; 'go-'; 'mo-' ;'yo-'; 'co-'];
figure;plot(1:nSeg,tp_w(4+2*nRad+1:4+2*nRad+nSeg),'ko-');hold on;
xlabel('Segment Index');ylabel('Activity Estimation');title('Constraint Effects');
strLegend_w=cell(nw+1,1);
strLegend_w{1}='truth';
for k=1:nw
    actL_w=pVals_w(k,4+2*nRad+1:4+2*nRad+nSeg);
    plot(1:nSeg,actL_w,formats_w(k,:));
    strLegend_w{k+1}=sprintf('Shape weight=%.1f',weights(k));   
end
legend(strLegend_w);
hold off;

figure; plot(weights,dscL_w,'b*-');
xlabel('weight');ylabel('DSC');title('Segmentation Results');

figure; plot(weights,fVals_w,'b*-');
xlabel('Weight');ylabel('fVal');title('Final Objective Function(Weight)');





