function [u,dx,dy,bx,by]=TV_SBM(u,Img,dx,dy,bx,by,Para)
[Nx,Ny]=size(u);
mu=Para.mu;lambda=Para.lambda;
%%update u
%center
A=dx(1:Nx-2,2:Ny-1)-dx(2:Nx-1,2:Ny-1)+dy(2:Nx-1,1:Ny-2)-dy(2:Nx-1,2:Ny-1);
A=-A+(bx(1:Nx-2,2:Ny-1)-bx(2:Nx-1,2:Ny-1)+by(2:Nx-1,1:Ny-2)-by(2:Nx-1,2:Ny-1));
B=u(3:Nx,2:Ny-1)+u(1:Nx-2,2:Ny-1)+u(2:Nx-1,3:Ny)+u(2:Nx-1,1:Ny-2);
G(2:Nx-1,2:Ny-1)=(lambda*(A+B)+mu*Img(2:Nx-1,2:Ny-1))./(4*lambda+mu);
%up border
A=dx(1,2:Ny-1)+dy(1,1:Ny-2)-dy(1,2:Ny-1);
A=-A+(bx(1,2:Ny-1)+by(1,1:Ny-2)-by(1,2:Ny-1));
B=u(2,2:Ny-1)+u(1,3:Ny)+u(1,1:Ny-2);
G(1,2:Ny-1)=(lambda*(A+B)+mu*Img(1,2:Ny-1))./(3*lambda+mu);
%bottom border
A=dx(Nx-1,2:Ny-1)-dx(Nx,2:Ny-1)+dy(Nx,1:Ny-2)-dy(Nx,2:Ny-1);
A=-A+(bx(Nx-1,2:Ny-1)-bx(Nx,2:Ny-1)+by(Nx,1:Ny-2)-by(Nx,2:Ny-1));
B=u(Nx-1,2:Ny-1)+u(Nx,3:Ny)+u(Nx,1:Ny-2);
G(Nx,2:Ny-1)=(lambda*(A+B)+mu*Img(Nx,2:Ny-1))./(3*lambda+mu);
%left border
A=dx(1:Nx-2,1)-dx(2:Nx-1,1)-dy(2:Nx-1,1);
A=-A+(bx(1:Nx-2,1)-bx(2:Nx-1,1)-by(2:Nx-1,1));
B=u(3:Nx,1)+u(1:Nx-2,1)+u(2:Nx-1,2);
G(2:Nx-1,1)=(lambda*(A+B)+mu*Img(2:Nx-1,1))./(3*lambda+mu);
%right border
A=dx(1:Nx-2,Ny)-dx(2:Nx-1,Ny)+dy(2:Nx-1,Ny-1)-dy(2:Nx-1,Ny);
A=-A+(bx(1:Nx-2,Ny)-bx(2:Nx-1,Ny)+by(2:Nx-1,Ny-1)-by(2:Nx-1,Ny));
B=u(3:Nx,Ny)+u(1:Nx-2,Ny)+u(2:Nx-1,Ny-1);
G(2:Nx-1,Ny)=(lambda*(A+B)+mu*Img(2:Nx-1,Ny))./(3*lambda+mu);
%cornors
G(1,1)=(mu*Img(1,1)+lambda*(u(2,1)+u(1,2)+dx(1,1)+dy(1,1)-bx(1,1)-by(1,1)))./(2*lambda+mu);
G(1,Ny)=(mu*Img(1,Ny)+lambda*(u(2,Ny)+u(1,Ny-1)-dx(1,Ny)+dy(1,Ny-1)-dy(1,1)+bx(1,Ny)-by(1,Ny-1)+by(1,Ny)))./(2*lambda+mu);
G(Nx,1)=(mu*Img(Nx,1)+lambda*(u(Nx-1,1)+u(Nx,2)-dx(Nx-1,1)+bx(Nx-1,1)-dx(Nx,1)+bx(Nx,1)-dy(Nx,1)+by(Nx,1)))./(2*lambda+mu);
G(Nx,Ny)=(mu*Img(Nx,Ny)+lambda*(u(Nx-1,1)+u(1,Ny-1)+dx(Nx-1,Ny)-dx(Nx,Ny)+dy(Nx,Ny-1)-dy(Nx,Ny)-(bx(Nx-1,Ny)...
    -bx(Nx,Ny)+by(Nx,Ny-1)-by(Nx,Ny))))./(2*lambda+mu);
u=G;

if Para.UseEdge==1
    gb=Para.Ge;
else
    gb=1;
end

uu = NeumannBoundCond(u); 
[ux,uy]=gradient(uu);
dx=ux+bx;dy=uy+by;dx=dx+1e-6*(dx==0);dy=dy+1e-6*(dy==0);
normd=sqrt(dx.^2+dy.^2);
dx=dx./normd.*max(normd-1./lambda.*gb,0);dy=dy./normd.*max(normd-1./lambda.*gb,0);
bx=bx+ux-dx;by=by+uy-dy;


function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]); 
