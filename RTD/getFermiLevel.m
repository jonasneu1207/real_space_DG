function [EfL, EfR] = getFermiLevel(MAT)
% getFermiLevel_nIntraband calculates the fermi level for conduction band
% quantum transport

% definition of the constants
h  = 6.626E-34;
hb = h/2/pi;
q  = 1.602E-19;
kB = 1.38E-23;
m0 = 9.109E-31;

dE  = 2E-4;

%% data decompression

T = MAT.Temp;
Vc = MAT.Vc;
Vh = MAT.Vh;
mc = MAT.mc;
mh = MAT.mlh;
N_D = MAT.N_d;
N_A = MAT.N_a;

beta = q/kB/T;

%% evalution of the Fermi-Level at the left contact

if N_D(1) > N_A(1)
    
    % energy range for the integration of the spectral carrier density
    EL   = min(Vc(1))-10*dE : dE : max(Vc(1))+0.3;
    % possible fermi level candidates
    mufL = min(Vc(1))-0.2 : 1E-3 : max(Vc(1))+0.5;
    
    % charge carriers left constact
    nL = getChargeCarrierDensity_nIntraband(EL, Vc(1), mufL, beta, mc(1));
    
    rho = nL/1E+6-N_D(1);
    [~, idx] = min(abs(log(rho)));
    
    EfL = mufL(idx);
    
else % N_D(1) < N_A(1)
    
    % energy range for the integration of the spectral carrier density
    EL   = min(Vh(1))-0.3 : dE : max(Vh(1))+10*dE;
    % possible fermi level candidates
    mufL = min(Vh(1))-0.5 : 1E-3 : max(Vh(1))+0.2;
    
    % charge carriers left constact
    pL = getChargeCarrierDensity_pIntraband(EL, Vh(1), mufL, beta, mh(1));
    
    rho = pL/1E+6-N_A(1);
    [~, idx] = min(abs(log(rho)));
    
    EfL = mufL(idx);
    
end

%% evalution of the Fermi-Level at the right contact


if N_D(end) > N_A(end)
    
    % energy range for the integration of the spectral carrier density
    ER   = min(Vc(end))-10*dE : dE : max(Vc(end))+0.3;
    % possible fermi level candidates
    mufR = min(Vc(end))-0.2 : 1E-3 : max(Vc(end))+0.5;
    
    % charge carriers left constact
    nR = getChargeCarrierDensity_nIntraband(ER, Vc(end), mufL, beta, mc(end));
    
    rho = nR/1E+6-N_D(end);
    [~, idx] = min(abs(log(rho)));
    
    EfR = mufR(idx);
    
else % N_D(end) < N_A(end)
    
    % energy range for the integration of the spectral carrier density
    ER   = min(Vh(end))-0.3 : dE : max(Vh(end))+10*dE;
    % possible fermi level candidates
    mufR = min(Vh(end))-0.5 : 1E-3 : max(Vh(end))+0.2;
    
    % charge carriers left constact
    pR = getChargeCarrierDensity_pIntraband(ER, Vh(end), mufR, beta, mh(end));
    
    rho = pR/1E+6-N_A(end);
    [~, idx] = min(abs(log(rho)));
    
    EfR = mufR(idx);
    
end