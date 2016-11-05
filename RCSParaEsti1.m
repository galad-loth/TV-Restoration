function [mu,alpha,var_n]=RCSParaEsti1(Img,W,L)
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
        Data2=SubImg.^2;Data4=SubImg.^4;Data6=SubImg.^6;
        m2=mean2(Data2);m4=mean2(Data4);m6=mean2(Data6);
%         nA=4*m2*(m4^2-m2*m6)/(2*m6*m2^2-m4*m6-m2*m4^2);
        nA=L;
        aA=(nA+1)*m2*m2/(nA*m4-(nA+1)*m2.^2);
        alpha(ii,jj)=aA;
        mu(ii,jj)=m2;
        
    end
end
alpha(alpha<=0)=10e3;
alpha(alpha>30)=30;
var_n=0;