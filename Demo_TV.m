clear all;close all;clc
%%prepare the image
Img=imread('Image/Roma.tif');
[Nx,Ny,Nc]=size(Img);
if Nc>1;Img=rgb2gray(Img);end
Img=double(Img(:,:,1));

Filter=fspecial('gaussian',5,1.2);
% ImgS=medfilt2(Img,[5 5]) ; 
ImgS=conv2(Img,Filter,'same');
figure(1);set(gcf,'position',[250 250 600 600]);
subplot(2,2,1);imagesc(Img);colormap(gray);
subplot(2,2,2);imagesc(ImgS);colormap(gray);

%% Parameter
u=Img; 
StopCret=0.2; u_old=0;
lam=0.01; NumIter=10; dt=0.2; eps=1; 

while (mean(mean(abs(u -u_old))) > StopCret),  % iterate until convergence
   u_old=u;
   u=tv(u,NumIter,dt,eps,lam,Img); % scalar lam
end % while

subplot(2,2,4);imagesc(u);colormap(gray);
