%%An example demo for GSSIP to generated neuron-shaped illumination
%%refernce "Highly efficient two-photon patterned photostimulation for cellular-resolution control of neuronal activity"
%%this demo perform that generated phase masks for speckle-free CGH form a
%%binarized target image

clc
clear
close all

cellNum = 1;    %select output target neuron No. (cellNum = 1 or 2 or 3)
load(['.\targetNeuronShape\neuron_target',num2str(cellNum),'.mat']);  % load binarized target pattern

%%optical paremeters 
N = 512;    %SLM pixels number
F = 200/20; %objective focal length
dxl = 15E-3;    %SLM pixels size
Magi_slm2foci = 1;  %Magnification factor between SLM and objective
liptimes = 20; %literations

lamda = 920e-6; %imput beam wavelength
ROT = 0;    %Angle between long axis (imput Gaussian beam) and horizontal direction
in_wx = 1.5;  %size of long axis imput Gaussian beam, mm
in_wy = 1.5;  %size of short axis imput Gaussian beam, mm

%%input beam sampling
in_M = 1;   %Scale factor
in_wmx = in_wx/dxl*in_M;
in_wmy = in_wy/dxl*in_M;
M = floor(5*2*max(in_wmx,in_wmy))+1;    %magnify 5 times to enable reasonal FFT
if mod(M,2) == 1
	M = M+1;
end

%%input and output setting
xs1 = -M/2:M/2-1; ys1=xs1; [XS1,YS1] = meshgrid(xs1,ys1);
XS1_RO = XS1.*cos((90-ROT)/180*pi)+YS1.*sin((90-ROT)/180*pi);
YS1_RO = YS1.*cos((90-ROT)/180*pi)-XS1.*sin((90-ROT)/180*pi);
U_in = exp(-(XS1_RO.^2./in_wmx^2+YS1_RO.^2./in_wmy^2));
% figure,imshow(U_in((M-N)/2+1:(M+N)/2,(M-N)/2+1:(M+N)/2),[]); colormap(jet);

imgSize = N/Magi_slm2foci*M/N;
target_Img = imresize(img,[imgSize,imgSize]);
% figure,imshow(target_Img,[]); colormap(jet);
U_out = zeros(M);  
U_out(floor((M-imgSize)/2)+1:floor((M+imgSize)/2),floor((M-imgSize)/2)+1:floor((M+imgSize)/2)) = target_Img;
Signal = imbinarize(U_out,0.3);
% figure,imshow(U_out,[]); colormap(jet);


%%GSSIP
img_size = N/Magi_slm2foci;
target_img = imresize(img,[img_size,img_size]);
u_out = zeros(N); %target img for calculation of inner and outer tangent circles
u_out(floor((N-img_size)/2)+1:floor((N+img_size)/2),floor((N-img_size)/2)+1:floor((N+img_size)/2)) = target_img;
u_outEdge = bwperim(imbinarize(img),8);
[x,y]  = find(u_outEdge==1);
edge(1,:) = x;
edge(2,:) = y;
[~,~,smallR,bigR] = getCircle(edge);
smallR = smallR/Magi_slm2foci*M/N;  %radius of inner tangent circles
bigR = bigR/Magi_slm2foci*M/N;  %radius of outer tangent circles
[t_GSSIP,targetALL_GSSIP,phaseALL_GSSIP] = GSSIP_Template(U_in, U_out, liptimes,in_wmx,in_wmy,smallR,bigR);   %gssip function
RMES_GSSIP = sqrt(sum(sum((targetALL_GSSIP(Signal==1)-U_out(Signal==1)).^2))/sum(sum((U_out(Signal==1)).^2))); %performence for gssip
yita_GSSIP = sum(sum(targetALL_GSSIP(Signal==1)))/sum(sum(targetALL_GSSIP));
target_GSSIP = targetALL_GSSIP(floor((M-N)/2)+1:floor((M+N)/2),floor((M-N)/2)+1:floor((M+N)/2));
phase_GSSIP = phaseALL_GSSIP(floor((M-N)/2)+1:floor((M+N)/2),floor((M-N)/2)+1:floor((M+N)/2));

%%GS
[t_GS,targetALL_GS,phaseALL_GS] = GS_Template(U_in, U_out, liptimes);
RMES_GS = sqrt(sum(sum((targetALL_GS(Signal==1)-U_out(Signal==1)).^2))/sum(sum((U_out(Signal==1)).^2))); %performence for gs
yita_GS = sum(sum(targetALL_GS(Signal==1)))/sum(sum(targetALL_GS));
target_GS = targetALL_GS(floor((M-N)/2)+1:floor((M+N)/2),floor((M-N)/2)+1:floor((M+N)/2));
phase_GS = phaseALL_GS(floor((M-N)/2)+1:floor((M+N)/2),floor((M-N)/2)+1:floor((M+N)/2));


%%slm phase mask 8bit RGB
center_x = N/2;center_y = center_x; %BNS center
locy = center_y-N/2+1 : center_y+N/2;
locx = center_x-N/2+1 : center_x+N/2;

dist_GSSIP = phase_GSSIP/2/pi*65536;
SLM_GSSIP = uint8(zeros(N,N,3)); %pre BNS-SLM 
SLM_GSSIP(locy, locx, 1) = uint8(mod(dist_GSSIP,256));    %R
SLM_GSSIP(locy, locx, 2) = uint8(floor(dist_GSSIP/256));  %G

dist_GS = phase_GS/2/pi*65536;
SLM_GS = uint8(zeros(N,N,3)); %pre BNS-SLM 
SLM_GS(locy, locx, 1) = uint8(mod(dist_GS,256));    %R
SLM_GS(locy, locx, 2) = uint8(floor(dist_GS/256));  %G

figure,
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.64]);
subplot(2,2,1),imshow(target_GSSIP,[]); colormap(hot); title({'Simulated Output Illumination by GSSIP';['RMSE = ', num2str(RMES_GSSIP),...
    'efficiency = ', num2str(yita_GSSIP)]});
subplot(2,2,2),imshow(target_GS,[]); colormap(hot); title({'Simulated Output Illumination by GS';['RMSE = ', num2str(RMES_GS),...
    'efficiency = ', num2str(yita_GS)]});
subplot(2,2,3),imshow(SLM_GSSIP,[]); title('phase mask by GSSIP');
subplot(2,2,4),imshow(SLM_GS,[]); title('phase mask by GSSIP');

