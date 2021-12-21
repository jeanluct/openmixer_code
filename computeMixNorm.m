function MixNorm = computeMixNorm(c, Nx, Ny)

% c must be of size Nx X Ny

% return Mc (size Nx X Ny) which is given as [M]c
fftc = abs(fft2(c)/(Nx*Ny)).^2;
Mc = zeros(Nx, Ny);

for kx = 0:(Nx/2-1)
for ky = 0:(Ny/2-1)
Mc(kx+1,ky+1)=fftc(kx+1,ky+1)/sqrt(1+4*pi^2*(kx^2+ky^2)); 

if (kx~=0)
Mc(Nx-kx+1,ky+1)=fftc(Nx-kx+1,ky+1)/sqrt(1+4*pi^2*(kx^2+ky^2));
end
 
if (ky~=0)
Mc(kx+1,Ny-ky+1)=fftc(kx+1,Ny-ky+1)/sqrt(1+4*pi^2*(kx^2+ky^2));
end

if ((kx~=0)&(ky~=0))
Mc(Nx-kx+1,Ny-ky+1)=fftc(Nx-kx+1,Ny-ky+1)/sqrt(1+4*pi^2*(kx^2+ky^2));
end

end

end

MixNorm = sum(sum(Mc));
