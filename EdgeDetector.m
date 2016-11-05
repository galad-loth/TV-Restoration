function Gb=EdgeDetector(Img,sigma_ED)
%%edge detector with gradient method
beta = 50/max(Img(:))^2;
[Ny,Nx]=size(Img);
Gaussian = fspecial('gaussian',[Ny Nx],sigma_ED);
ImgS = ifftshift( ifft2( fft2(Img) .* fft2(Gaussian) ) );
% Gaussian = fspecial('gaussian',3,sigma_ED);
% ImgS=conv2(Img,Gaussian,'same');
% Norm of the gradient of the smoothed image
[Gx,Gy] = gradient(ImgS);
NormGrad = sqrt(Gx.^2 + Gy.^2);
% Edge detector function
Gb = ones(Ny,Nx)./ (1 + beta* NormGrad.^2);

