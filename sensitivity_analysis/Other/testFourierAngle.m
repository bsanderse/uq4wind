%%
clear all
close all

%% generate some sines

% uncertain parameter
z1 = 0.2;
z2 = 0.9;
Nz  = 8;

z = z1 + (z2-z1)*rand(Nz,1);

x1 = 0;
x2 = 1;
Nx = 30;
x  = linspace(x1,x2,Nx);

Fn = zeros(Nx,Nz);

Famp   = zeros(Nz,1);
Fphase = zeros(Nz,1);

Freal  = zeros(Nz,1);
Fimag  = zeros(Nz,1);

%%
for k=1:Nz
    
    % create physical signal
    Fn(:,k) = 0.1 + z(k).*cos(2*pi*x - 4*pi*z(k).^2 + 3);
    
    % get Fourier transform
    [Fhat,Fabs,Fphi]  = getFourierCoefficientsRealData(Fn,2);
    
    % get amplitude and phase of second mode
    Famp(k) = 2*abs(Fhat(2,k));
    Fphase(k) = angle(Fhat(2,k));
    
    % or use the real and imaginary coefficients
    Freal(k) = real(Fhat(2,k));
    Fimag(k) = imag(Fhat(2,k));
    
end

%%
figure
plot(Fn)

[z,sort_id] =sort(z);

%% amplitude
figure
plot(z,Famp(sort_id),'o-')
hold on
plot(z,z,'d-');
title('amplitude')
legend('amplitude (derived)','amplitude (given)')

%% phase
figure
plot(z,Fphase(sort_id),'o-');
hold on
plot(z,unwrap(Fphase(sort_id)),'s-');
% exact angle:
plot(z,-4*pi*z.^2+3,'d');
legend('phi (derived)','phi (derived, unwrapped)','phi (given)')

%% real and imaginary components
figure
plot(z,Freal(sort_id),'o-')
hold on
plot(z,0.5*z.*cos(-4*pi*z.^2+3),'d-')

plot(z,Fimag(sort_id),'s-')
% abs(z)*exp(i*theta) = abs(z)*(cos(theta) + i*sin(theta))
% = a + b*i
% a = abs(z)*cos(theta)
% b = abs(z)*sin(theta)
plot(z,0.5*z.*sin(-4*pi*z.^2+3),'d-')
title('real and imaginary components')
legend('real component (derived)','real component (given)','imaginary component (derived)','imaginary component (given)');
