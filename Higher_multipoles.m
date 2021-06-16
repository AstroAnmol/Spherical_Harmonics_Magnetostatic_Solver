clear all;
%% Inputs

B0 = 6e-4;
susc = 0.96; % Magnetic susceptibility
a = 1.4e-6;  % Grain radius, meters
sep=2.2;
alpha=0;
L=10;       %Number of multipoles used

%% Code to solve for coefficients

mu0 = 4*pi*1e-07;
mu = (1+susc)*mu0;
Hmag = B0/(mu0);  % Applied magnetic field, A/m
sep=sep*a;

R1=[0 0 0]';
R2=[0 0 sep]';

alpha = alpha*pi/180; % change angle into radians
H0 = [Hmag*sin(alpha) 0 Hmag*cos(alpha)]'; % A/m

% parallel and perpendicular components of the applied magnetic field
H_perp=H0(1);
H_prll=H0(3);

% Creating the L X L matrices 
X=zeros(L,L);
Delta0=zeros(L,L);
Gamma0=zeros(L,L);
Delta1=zeros(L,L);
Gamma1=zeros(L,L);
for i=1:L
    for j=1:L
        % X matrix
        if i==j
            X(i,j)= i*(mu/mu0) + i +1;
        end
        % Delta and Gamma matrix for m=0
        Delta0(i,j)= ((-1)^(i))*(i*(mu/mu0)-i)*nchoosek(i+j,j)*...
            (a^(2*i+1))/(sep^(i+j+1));
        Gamma0(i,j)= ((-1)^(i+j))*Delta0(i,j);
        
        % Delta and Gamma matrix for m=1
        Delta1(i,j)= ((-1)^(i+1))*(i*(mu/mu0)-i)*nchoosek(i+j,j+1)*...
            (a^(2*i+1))/(sep^(i+j+1));
        Gamma1(i,j)= ((-1)^(i+j))*Delta1(i,j);
    end
end

% q vectors for m =0,1
q0=zeros(L,1);
q1=zeros(L,1);

q0(1)=H_prll*(a^3)*(1-mu/mu0);
q1(1)=-H_perp*(a^3)*(1-mu/mu0);

% creating the 2L X 2L matrix

A0=zeros(2*L);
A1=zeros(2*L);

A0(1:L,1:L)=X;
A0(L+1:end,1:L)=Gamma0;
A0(1:L,L+1:end)=Delta0;
A0(L+1:end,L+1:end)=X;


A1(1:L,1:L)=X;
A1(L+1:end,1:L)=Gamma1;
A1(1:L,L+1:end)=Delta1;
A1(L+1:end,L+1:end)=X;

Q0=[q0;q0];
Q1=[q1;q1];

Beta0=A0\Q0;
Beta1=A1\Q1;

Beta1_0=Beta0(1:L);
Beta2_0=Beta0(L+1:2*L);
Beta1_1=Beta1(1:L);
Beta2_1=Beta1(L+1:2*L);

%% Computing the Magnetic Field

%Legendre Polynomial for different values
n=500;
dx = pi/60;
col = 0:dx:2*pi;
r1 = linspace(a, 16*a, n);
[theta,R] = meshgrid(col,r1);


Hrs=0;
Hr=0;

for l=1:L
    for m=0:1
        for s=m:L
            Psm=legendre(s,cos(theta));
            if s~=0
                Psm=reshape(Psm(m+1,:,:),size(R));
            end
            Hrs=Hrs+ (-1)^(s+m) * nchoosek(l+s,s+m)*s.*(R.^(s-1)).*Psm/...
                (sep^(l+s+1));
        end
        if m==0
            Plm=legendre(l,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Hr=Hr+(l+1)*Beta1_0(l).*Plm./(R.^(l+2)) - Beta2_0(l)*Hrs;
        elseif m==1
            Plm=legendre(l,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Hr=Hr+(l+1)*Beta1_1(l).*Plm./(R.^(l+2)) - Beta2_1(l)*Hrs;
        end
    end
end


Hths=0;
Hth=0;

for l=1:L
    for m=0:1
        for s=m:L
            Psm=legendre(s,cos(theta));
            Ps1m=legendre(s+1,cos(theta));
            if s~=0
                Psm=reshape(Psm(m+1,:,:),size(R));
            end
            Ps1m=reshape(Ps1m(m+1,:,:),size(R));
            %compute derivative of the associated legendre function
            dPsm=((m-s-1).*Ps1m + (s+1).*cos(theta).*Psm)./(-sin(theta));
            Hths=Hths+ (-1)^(s+m) * nchoosek(l+s,s+m).*(R.^(s-1)).*dPsm/...
                (sep^(l+s+1));
        end
        if m==0
            Plm=legendre(l,cos(theta));
            Pl1m=legendre(l+1,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Pl1m=reshape(Pl1m(m+1,:,:),size(R));
            
            %compute derivative of the associated legendre function
            dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
            Hth=Hth+Beta1_0(l).*Plm./(R.^(l+2)) - Beta2_0(l)*Hths;
        elseif m==1
            Plm=legendre(l,cos(theta));
            Pl1m=legendre(l+1,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Pl1m=reshape(Pl1m(m+1,:,:),size(R));
            
            %compute derivative of the associated legendre function
            dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
            Hth=Hth+Beta1_1(l).*Plm./(R.^(l+2)) - Beta2_1(l)*Hths;
        end
    end
end

%% Formulating the Maxwell Stress Tensor in Spherical Coordinates
size=size(Hr);
for i=1:size(1)
    for j=1:size(2)
        H=[Hr(i,j);Hth(i,j)];
        h=norm(H);
        sigma=(mu0.*H*H' - 0.5*(h^2).*eye(2));
    end
end
