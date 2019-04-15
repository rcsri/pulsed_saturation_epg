function g = lineshape_superlorentzian(f, T2)

u = linspace(0,1,128);
c = abs(3*u.^2 - 1);
g = sqrt(2/pi) * T2 * trapz(u, exp(-2*(2*pi*f*T2./c).^2) ./ c); % Super-Lorentzian (Morrison et al 1995 Eq 8)