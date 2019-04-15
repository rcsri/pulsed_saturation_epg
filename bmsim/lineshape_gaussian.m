function g = lineshape_gaussian(f, T2)

g = (T2 / sqrt(2*pi)) * exp ( (-(2*pi*f*T2).^2) / 2); % Gaussian (Morrison et al 1995 Eq. 5)
