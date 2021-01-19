function R=rotxn(flip,phi,w1,dw0);
w=sqrt(w1*w1+dw0*dw0);
nx=w1/w*cos(phi);
ny=w1/w*sin(phi);
nz=dw0/w;
t=flip/w1;
theta=t*w;
R11=nx*nx*(1-cos(theta))+cos(theta);
R12=nx*ny*(1-cos(theta))+nz*sin(theta);
R13=nx*nz*(1-cos(theta))-ny*sin(theta);
R21=nx*ny*(1-cos(theta))-nz*sin(theta);
R22=ny*ny*(1-cos(theta))+cos(theta);
R23=ny*nz*(1-cos(theta))+nx*sin(theta);
R31=nx*nz*(1-cos(theta))+ny*sin(theta);
R32=ny*nz*(1-cos(theta))-nx*sin(theta);
R33=nz*nz*(1-cos(theta))+cos(theta);
R=[R11 R12 R13;R21 R22 R23;R31 R32 R33];

