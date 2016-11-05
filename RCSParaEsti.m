function [mu,alpha,var_n]=RCSParaEsti(Img,W,L)
[Nx,Ny]=size(Img);
alpha=zeros(Nx,Ny);mu=alpha;var_n=alpha;
for ii=1:Nx
    for jj=1:Ny
        if ii<(W+1)/2
            idxx=1:ii+(W-1)/2;
        elseif ii>Nx-(W-1)/2
            idxx=ii-(W-1)/2:Nx;
        else
            idxx=ii-(W-1)/2:ii+(W-1)/2;
        end
        if jj<(W+1)/2
            idxy=1:jj+(W-1)/2;
        elseif jj>Ny-(W-1)/2
            idxy=jj-(W-1)/2:Ny;
        else
            idxy=jj-(W-1)/2:jj+(W-1)/2;
        end
        SubImg=Img(idxx,idxy);
        mu_temp=sum(sum(SubImg))/W/W;
        var=sum(sum((SubImg-mu_temp).^2))/W/W;
        var_n(ii,jj)=var/mu_temp/mu_temp;
        alpha(ii,jj)=(1+1/L)/(var_n(ii,jj)-1/L);
        mu(ii,jj)=mu_temp;
    end
end