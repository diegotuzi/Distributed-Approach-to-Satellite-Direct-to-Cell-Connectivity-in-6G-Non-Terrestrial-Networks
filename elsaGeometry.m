%% arrayGeometry_lsaSt
% Diego Tuzi UniBwM (2022)

% N: number of elements
% fc: frequency
% d: one-dimensional linear spacing in wavelenghts
% stCoeff: spatial tapering coefficient
% stCoeff= 0.1 optimized for fc=2Ghz and d=6.67lambda 
% stCoeff=1 for LSA as in US patent 6433754

function [x,y]=arrayGeometry_lsaSt(N,d,stCoeff)
    n=1:N;
    dn=(d*(1-stCoeff)/N)*n+d*stCoeff; % modified one-dim. linear spacing in w.
    
    r=dn.*sqrt(n/pi); % modified polar equation - radius
    
    tau=(1+sqrt(5))/2; % golden angle
    
    theta=2*pi*n*tau; % polar equation - angle
    
    % conversion to cartesian coordinates
    x=r.*cos(theta);
    y=r.*sin(theta);
end