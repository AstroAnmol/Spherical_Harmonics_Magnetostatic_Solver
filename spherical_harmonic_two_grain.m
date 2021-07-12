function f=spherical_harmonic_two_grain(B0, susc, a, sep, alpha, L)

debug=0;

mu0 = 4*pi*1e-07;
mu = (1+susc)*mu0;
Hmag = B0/(mu0);  % Applied magnetic field, A/m
sep=sep*a;

R1=[0 0 sep]';
R2=[0 0 0]';

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

% %Using SVD to solve the linear system
% [U0,S0,V0] = svd(A0);
% A0_inv=V0*inv(S0)*U0';
% Beta0_dash=A0_inv*Q0;


Beta0=A0\Q0;
Beta1=A1\Q1;

Beta1_0=Beta0(1:L);
Beta2_0=Beta0(L+1:2*L);
Beta1_1=Beta1(1:L);
Beta2_1=Beta1(L+1:2*L);

%% Computing the Magnetic Field

%Legendre Polynomial for different values
dang = pi/180;
col = dang:dang:2*pi;
az = dang:dang:pi;
dr=a/100;
r1=a-dr:dr:a+20*dr;
[phi,theta,R] = meshgrid(az,col,r1);


Hrs=0;
Hr=0;
Hths=0;
Hth=0;
Hphis=0;
Hphi=0;

for l=1:L
    for m=0:1
        Hrs=0;
        Hths=0;
        
        for s=m:L
            Psm=legendre(s,cos(theta));
            Ps1m=legendre(s+1,cos(theta));
            if s~=0
                Psm=reshape(Psm(m+1,:,:),size(R));
            end
            Ps1m=reshape(Ps1m(m+1,:,:),size(R));
            %compute derivative of the associated legendre function
            dPsm=((m-s-1).*Ps1m + (s+1).*cos(theta).*Psm)./(-sin(theta));
            % R component
            Hrs=Hrs+ (-1)^(s+m) * nchoosek(l+s,s+m)*s.*(R.^(s-1)).*Psm/...
                (sep^(l+s+1));
            % Theta component
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
            % Theta Component
            Hth=Hth+...
                (Beta1_0(l).*dPlm./(R.^(l+2)) + Beta2_0(l)*Hths).*cos(m*phi);
            % R Component
            Hr=Hr+...
                ((l+1)*Beta1_0(l).*Plm./(R.^(l+2)) -...
                Beta2_0(l)*Hrs).*cos(m*phi);
        elseif m==1
            Plm=legendre(l,cos(theta));
            Pl1m=legendre(l+1,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Pl1m=reshape(Pl1m(m+1,:,:),size(R));
            
            %compute derivative of the associated legendre function
            dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
            % Theta Component
            Hth=Hth+...
                (Beta1_1(l).*dPlm./(R.^(l+2)) + Beta2_1(l)*Hths).*cos(m*phi);
            % R Component
            Hr=Hr+...
                ((l+1)*Beta1_1(l).*Plm./(R.^(l+2)) -...
                Beta2_1(l)*Hrs).*cos(m*phi);
        end
    end
end
%Phi component
for l=1:L
    Hphis=0;
    for s=1:L
        Ps1=legendre(s,cos(theta));
        Ps1=reshape(Ps1(2,:,:),size(R));
        Hphis=Hphis+ (-1)^(1+s) * nchoosek(l+s,s+1).*(R.^(s-1)).*Ps1...
            ./ (sin(theta).*sep^(l+s+1));
    end
    Pl1=legendre(l,cos(theta));
    Pl1=reshape(Pl1(2,:,:),size(R));
    Hphi=Hphi+...
                ((l+1)*Beta1_1(l).*Pl1./(R.^(l+2)) -...
                Beta2_1(l)*Hphis).*sin(phi);
end
%% plot Magnetic Field
if(debug==1)
%[x,y,z]=sph2cart(phi,theta, R);
x=R.*sin(theta).*sin(phi);
y=R.*cos(phi);
z=R.*cos(theta).*sin(phi);
% pre=[sin(ph)*cos(th), cos(ph)*cos(th), -sin(th);...
%     sin(ph)*sin(th), cos(th)*sin(th),  cos(th);...
%     cos(ph),         -sin(ph),         0];

Hmat = sqrt(Hr.^2+Hth.^2+Hphi.^2);
figure;
Ang=1*360/4;
pc=pcolor(squeeze(z(:,Ang,:))./a,squeeze(x(:,Ang,:))./a,squeeze(Hmat(:,Ang,:))); set(pc, 'EdgeColor', 'none');
xlim([-3 3]); ylim([-3 3]);
title('|H|');
xlabel('z');
ylabel('x');
colorbar;
axis equal;
grid on;
xlim([-3 3]); ylim([-3 3]);
end
%% Formulating the Maxwell Stress Tensor in Spherical Coordinates
Size_col=size(col);
Size_az=size(az);
fx=0;
fy=0;
fz=0;
for i=1:Size_col(2)
    if i==1 || i==Size_col(2)
        p=1;
    elseif mod(i,2)~=0
        p=2;
    else
        p=4;
    end
    for j=1:Size_az(2)
        if j==1 || j==Size_az(2)
            q=1;
        elseif mod(j,2)~=0
            q=2;
        else
            q=4;
        end
            th=col(i);
            ph=az(j);
            %transformation matrix
                pre=[sin(ph)*cos(th), cos(ph)*cos(th), -sin(th);...
                    sin(ph)*sin(th), cos(th)*sin(th),  cos(th);...
                    cos(ph),         -sin(ph),         0];
            H_sph=[Hr(i,j,2);Hphi(i,j,2);Hth(i,j,2)];% + CART2SPH*H0;
            H_cart=pre*H_sph;
            h=norm(H_sph);
            %Hmat(i,j)=h;
            Txr=mu0*H_cart(2)*H_sph(1);
            Tyr=mu0*H_cart(3)*H_sph(1);
            Tzr=mu0*H_cart(1)*H_sph(1);
            fx=fx+a*a*Txr*sin(ph)*p*q;
            fy=fy+a*a*Tyr*sin(ph)*p*q;
            fz=fz+a*a*Tzr*sin(ph)*p*q;
            
    end
end
f=[fx, fy, fz];
f=f*dang*dang/9;
end