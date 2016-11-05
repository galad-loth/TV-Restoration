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
g=Img;
mu=0.03;lambda=1;tol=0.01;
u=SplitBregmanROF(g,mu,lambda,tol);
subplot(2,2,4);imagesc(u);colormap(gray);
