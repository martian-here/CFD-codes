clc 
clear all 
close all

% Input
u_0 = 1;
imax = 51; jmax = 51; 
Re = 400; mu = 1/Re;
omega(1:imax,1:jmax) = 0;
si(1:imax,1:jmax) = 0;
u(1:imax,1:jmax) = 0; v(1:imax,1:jmax) = 0;
dx = 1/(imax-1); dy = 1/(jmax-1);
error(4,imax,jmax) = 0; epsilon = 0.001;
% IC & BC
u(:,jmax) = 1;

%
for k = 1:5000
    for i = 2:imax-1
        for j = 2:jmax-1
            u_old = u(i,j);
            v_old = v(i,j);
            u(i,j) = (si(i,j+1)-si(i,j-1))/(2*dy);
            v(i,j) = -(si(i+1,j)-si(i-1,j))/(2*dx);
            error(1,i,j) = abs(u(i,j)-u_old);
            error(2,i,j) = abs(v(i,j)-v_old);
            si_old = si(i,j);
            si(i,j) = (si(i+1,j)+si(i-1,j)+si(i,j+1)+si(i,j-1)+omega(i,j)*dx^2)/4;
            error(3,i,j) = abs(si(i,j)-si_old);
            omega_old = omega(i,j);
            ae = -min(u(i,j),0)+mu/dx; aw = max(u(i,j),0)+mu/dx;
            an = -min(v(i,j),0)+mu/dx; as = max(v(i,j),0)+mu/dx; 
            ap = ae+aw+an+as;
            omega(i,j) = (ae*omega(i+1,j) + aw*omega(i-1,j) + an*omega(i,j+1) + as*omega(i,j-1))/ap;
            error(4,i,j) = abs(omega(i,j)-omega_old);
        end
    end
  for i = 2:imax-1
        omega_old = omega(1,i);
        omega(1,i) = -2*si(2,i)/dx^2;
        error(4,1,i) = abs(omega(1,i)-omega_old);
        omega_old = omega(imax,i);
        omega(imax,i) = -2*si(imax-1,i)/dx^2;
        error(4,imax,i) = abs(omega(imax,i)-omega_old);
        omega_old = omega(i,1);
        omega(i,1) = -2*si(i,2)/dy^2;
        error(4,i,1) = abs(omega(i,1)-omega_old);
        omega_old = omega(i,imax);
        omega(i,imax) = -(2*si(i,imax-1)+2*u_0*dy)/dy^2;
        error(4,i,imax) = abs(omega(i,imax)-omega_old);
    end
    if max(error) < epsilon
        break
    end
end

k
xv = [1.0000 0.0000; 0.9688 -0.05906;0.9609 -0.07391;
    0.9531 -0.08864; 
    0.9453  -0.10313;   
 0.9063 -0.16914;
 0.8594 -0.22445;
 0.8047 -0.24533;   
 0.5 0.05454;    
 0.2344 0.17527;    
 0.2266 0.17507;    
 0.1563 0.16077;    
 0.0938 0.12317;    
 0.0781 0.1089;     
 0.0703 0.10091;    
 0.0625 0.09233;    
 0.0000 0.00000;];

yu = [1        1.0000;    
 0.9766   0.84123;   
 0.9688   0.78871;   
 0.9609   0.73722;   
 0.9531   0.68717;   
 0.8516   0.23151;   
 0.7344   0.00332;   
 0.6172  -0.136641;  
 0.5     -0.20581;  
 0.4531  -0.2109;   
 0.2813  -0.15662;  
 0.1719  -0.1015;   
 0.1016  -0.063434;  
 0.0703  -0.04775;  
 0.0625  -0.04192;  
 0.0547  -0.03717;  
 0.0000   0.00000; ];


xv_(1:17) = xv(17:-1:1,2);
yu_(1:17) = yu(17:-1:1,2);

x(1:17) = xv(17:-1:1,1);
scatter(x,xv_(:))
hold on
y = linspace(0,1,imax);
plot(y,v(:,(imax-1)/2+1))
legend('ghia et. el.','51*51 grid')
xlabel('X')
ylabel('v')
figure

x = yu(17:-1:1,1);
scatter(yu_(:),x)
hold on
y = linspace(0,1,imax);
plot(u((imax-1)/2+1,:),y)
xlabel('u')
ylabel('Y')
legend('ghia et. el.','51*51 grid')
figure

%
for i = 1:imax
    for j = 1:jmax
        p_omega(j,i) = omega(i,j);
        p_si(j,i) = si(i,j);
        p_u(j,i) = u(i,j);
        p_v(j,i) = v(i,j);
    end
end

contourf(p_si(:,:))
colormap("jet")
colorbar
hold on
quiver(p_u,p_v)
figure
contourf(p_omega(:,:))
colormap("jet")
colorbar
legend ('vorticity')
figure
a = linspace(0,1,51);
b = linspace(0,1,51);
quiver(a,b,p_u,p_v)
xlim([0 1])
ylim([0 1])


