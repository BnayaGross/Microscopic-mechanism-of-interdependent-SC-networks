function duout = fft2d(t,u)

v = size(u');
Nx = v(1);
Ny = v(2);
uhat = fft2(u);

  if rem(Nx, 2) ==0

  kap_x = (2*pi/Nx)*[-Nx/2:Nx/2-1];
   kap1x = fftshift(kap_x)';
 else
  kap_x = (2*pi/Nx)*[-(Nx-1)/2:(Nx-1)/2]; 
  kap1x = fftshift(kap_x);
  kap1x=[0 kap1x(1:end-1)]';
  end

  if rem(Ny, 2) ==0
  kap_y = (2*pi/Ny)*[-Ny/2:Ny/2-1];
  kap1y = fftshift(kap_y)';
else
  kap_y = (2*pi/Ny)*[-(Ny-1)/2:(Ny-1)/2]; 
  kap1y = fftshift(kap_y);
  kap1y=[0 kap1y(1:end-1)]';
 end
  

a1 = exp((2*t*(cos(kap1x)-1)))*exp((2*t*(cos(kap1y)-1)))';


duhat = a1.*uhat';



duout =ifft2(duhat');