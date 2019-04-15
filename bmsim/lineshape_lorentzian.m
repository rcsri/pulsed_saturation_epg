function g = lineshape_lorentzian(f, T2)

g = (T2 / pi) * (1/(1 + (2*pi*f*T2).^2)); % Lorentzian (Morrison et al 1995 Eq. 4)
