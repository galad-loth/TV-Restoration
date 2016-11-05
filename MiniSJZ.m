function z=MiniSJZ(g,u,d,mu,L)
[Nx,Ny]=size(g);
z=zeros(Nx,Ny);
z_old=ones(Nx,Ny);
while max(max(abs(z-z_old)))>0.01
    z_old=z;
    fz=z+exp(g-z)+mu/2/L*(z-u-d).^2;
    fdz=1-exp(g-z)+mu/L*(z-u-d);
    fddz=exp(g-z)+mu/L;
    z=z-fdz./(fddz+eps);
end