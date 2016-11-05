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
L=32; 
[muInt,nu,var_n]=RCSParaEsti1(Img,5,L); 
g=log(Img.^2+eps);
subplot(2,2,3);imagesc(g);colormap(gray)
mu=0.03;lambda=.5;tol=0.01;mu1=15;
u=g;d=zeros(Nx,Ny);
z=zeros(Nx,Ny);z_old=ones(Nx,Ny);
counter=0;
while sum(sum(abs(z-z_old)))>0.01
    counter=counter+1
    z_old=z;
    z=MiniSJZ1(g,u,d,mu1,L,nu,muInt);
    u_temp=z-d;
    u_temp_max=max(max(u_temp));
    u_temp=u_temp*255/u_temp_max;
    u=SplitBregmanROF(u_temp,mu,lambda,tol);
    subplot(2,2,3);imagesc(u);colormap(gray)
    u=u*u_temp_max/255;
    d=d+u-z;
end
subplot(2,2,4);imagesc(sqrt(exp(z)));colormap(gray)

