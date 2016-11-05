clear all;close all;clc
%%prepare the image
Img=imread('Image/Roma.tif');
[Nx,Ny,Nc]=size(Img);
if Nc>1;Img=rgb2gray(Img);end
Img=double(Img(:,:,1));

Filter=fspecial('gaussian',5,1.2);
% ImgS=medfilt2(Img,[5 5]) ; 
ImgS=conv2(Img,Filter,'same');
figure(1);set(gcf,'position',[250 150 500 500]);
subplot(2,2,1);imagesc(Img);colormap(gray);
subplot(2,2,2);imagesc(ImgS);colormap(gray);

%%
% g=log(Img+eps);
% subplot(2,2,3);imagesc(g);colormap(gray)
u=Img; 
StopCret=0.2; u_old=0;
lam=15; NumIter=50; dt=0.2; eps=1;L=32; counter=0;
[muInt,nu,var_n]=RCSParaEsti1(Img,5,L);  
while (mean(mean(abs(u -u_old))) > StopCret),  % iterate until convergence
   counter=counter+1
   u_old=u;
   u=tv_speckle_CorNeibor(u,NumIter,dt,eps,lam,nu,L,Img,5); % scalar lam
end % while
subplot(2,2,4);imagesc(u);colormap(gray);

