clear all
close all
clc

%% Dados da Aeronave

ftm = 0.3048;
lbkg = 0.4536;
slugkg = 14.5939;

S = 174*ftm^2;
b = 36*ftm;

g = 9.8065;
m = 2650*lbkg;
rho = 1.0555;
u0 = 220.1*ftm;
q = 0.5*(rho*(u0^2));       %Pressao Dinamica
cmedia = 4.9*ftm;

Ixx = 948*slugkg*ftm^2;
Iyy = 1346*slugkg*ftm^2;
Izz = 1967*slugkg*ftm^2;
Ixz = 0*slugkg*ftm^2;

ix = Ixz/Ixx;
iz = Ixz/Izz;

Cd0 = 0.0270;
Cdu = 0;
Cda = 0.121;
Ctxu = -0.096;
Cl0 = 0.307;
Clmax = 2.32; % not ok
Clu = 0;
Cla = 4.41;
Claponto = 1.7;
Clq = 3.9;
Cm0 = 0.04;
Cmu = 0;
Cma = -0.613;
Cmaponto = -7.27;
Cmq = -12.4;
Cmtu = 0;
Cmta = 0;

Cdde = 0;
Clde = 0.43;
Cmde = -1.122;
Cha = -0.0584;
Chde = -0.585;

Clb = -0.0923;
Clp = -0.484;
Clr = 0.0798;
Cyb = -0.393;
Cyp = -0.075;
Cyr = 0.214;
Cnb = 0.0587;
Cntb = 0;
Cnr = -0.0278;
Cnp = -0.0937;
% ok

Clda = 0.229;
Cldr = 0.0147;
Cyda = 0;
Cydr = 0.187;
Cnda = -0.0216;
Cndr = -0.0645;

%% Derivadas dimensionais de estabilidade Longitudinal

Xu = (-1)*(q*S*(Cdu + 2*Cd0))/(m*u0);
Xtu = (q*S*Ctxu)/(m*u0);
Xa = (-1)*(q*S*(Cda - Cl0))/m;
Zu = (-1)*(q*S*(Clu + 2*Cl0))/(m*u0);
Za = (-1)*(q*S*(Cla + 2*Cd0))/m;
Zaponto = (-1)*(q*S*cmedia*Claponto)/(2*m*u0);
Zq = (-1)*(q*S*cmedia*Clq)/(2*m*u0);
Mu = (q*S*cmedia*Cmu)/(Iyy*u0);
Mtu = (q*S*cmedia*Cmtu)/(Iyy*u0);
Ma = (q*S*cmedia*Cma)/(Iyy);
Mta = (q*S*cmedia*Cmta)/Iyy;
Maponto = (q*S*cmedia*cmedia*Cmaponto)/(2*Iyy*u0);
Mq = (q*S*cmedia*cmedia*Cmq)/(2*Iyy*u0);

Xde = (-1)*(q*S*Cdde)/m;
Zde = (-1)*(q*S*Clde)/m;
Mde = (-1)*(q*S*cmedia*Cmde)/Iyy;

%% Matrizes Longitudinais

Along(1,1) = Xu + Xtu;
Along(1,2) = Xa/u0;
Along(1,3) = 0;
Along(1,4) = -g;
Along(2,1) = (u0*Zu)/(u0 - Zaponto);
Along(2,2) = Za/(u0-Zaponto);
Along(2,3) = (u0*(u0 + Zq))/(u0 - Zaponto);
Along(2,4) = 0;
Along(3,1) = Mu + Mtu;
Along(3,2) = ((Ma + Mta)/u0) + ((Maponto * Za)/(u0*(u0 - Zaponto)));
Along(3,3) = ((Maponto*(Zq + u0))/(u0 - Zaponto)) + Mq;
Along(3,4) = 0;
Along(4,1) = 0;
Along(4,2) = 0;
Along(4,3) = 1;
Along(4,4) = 0;

Blong(1,1) = Xde;
Blong(2,1) = Zde/(u0 - Zaponto);
Blong(3,1) = Mde + ((Maponto*Zde)/(u0 - Zaponto));
Blong(4,1) = 0;

%% Matrizes C e D Longitudinais
Clong = eye(4);
Dlong = zeros (4,1);

%% Derivadas dimensionais de estabilidade Latero-Direcional

Yb = q*S*Cyb/m;
Yp = q*S*b*Cyp/(2*m*u0);
Yr = q*S*b*Cyr/(2*m*u0);
Lb = q*S*b*Clb/Ixx;
Lp = q*S*b*b*Clp/(2*Ixx*u0);
Lr = q*S*b*b*Clr/(2*Ixx*u0);
Nb = q*S*b*Cnb/Izz;
Ntb = q*S*b*Cntb/Izz;
Np = q*S*b*b*Cnp/(2*Izz*u0);
Nr = q*S*b*b*Cnr/(2*Izz*u0);

Yda = q*S*Cyda/m;
Ydr = q*S*Cydr/m;
Lda = q*S*Clda/Ixx;
Ldr = q*S*Cldr/Ixx;
Nda = q*S*Cnda/Izz;
Ndr = q*S*Cndr/Izz;

%% Matrizes Latero-Direcionais

Ald(1,1) = Yb/u0;
Ald(1,2) = Yp;
Ald(1,3) = Yr - u0;
Ald(1,4) = g;
Ald(2,1) = (Lb + ix*Nb +ix*Ntb)/(u0*(1- ix*iz));
Ald(2,2) = (Lp + ix*Np)/(1 - ix*iz);
Ald(2,3) = (Lr + ix*Nr)/(1 - ix*iz);
Ald(2,4) = 0;
Ald(3,1) = (iz*Lb + Nb + Ntb)/(u0*(1 - ix*iz));
Ald(3,2) = (iz*Lp + Np)/(1 - ix*iz);
Ald(3,3) = (iz*Lr + Nr)/(1 - ix*iz);
Ald(3,4) = 0;
Ald(4,1) = 0;
Ald(4,2) = 1;
Ald(4,3) = 0;
Ald(4,4) = 0;

Bld(1,1) = Yda;
Bld(1,2) = Ydr;
Bld(2,1) = (Lda + ix*Nda)/(1-ix*iz);
Bld(2,2) = (Ldr + ix*Ndr)/(1 - ix*iz);
Bld(3,1) = (iz*Lda + Nda)/(1 - ix*iz);
Bld(3,2) = (iz*Ldr + Ndr)/(1 - ix*iz);
Bld(4,1) = 0;
Bld(4,2) = 0;


%% Matrizes C e D L�tero-Direcionais
Cld = eye(4);
Dld = zeros (4,2);

%% Determina��o dos sistemas
syslong = ss(Along,Blong,Clong,Dlong);


sysld = ss(Ald,Bld,Cld,Dld);