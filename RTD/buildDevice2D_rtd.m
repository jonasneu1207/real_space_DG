function mat = buildDevice2D_rtd(device, Ly, dy)

switch device
    
    case 'n-RTD-constM-comp'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/2;

        Lc1 = ceil(19/a)*a;
        Ls1 = ceil(6/a)*a;
        Lb1 = ceil(4/a)*a;
        Lw  = ceil(6/a)*a;
        Lb2 = ceil(4/a)*a;
        Ls2 = ceil(6/a)*a;
        Lc2 = ceil(19/a)*a; 
        Ndb = 2E+18;
        Ndw = 0*1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'n-RTD-constM-comp-UDS-s1'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = 2*a;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 
        
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'n-RTD-constM-comp-UDS-s2'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 

        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'n-RTD-constM-comp-UDS-s3'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/2;
        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 

        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

            case 'n-RTD-constM-comp-UDS-s4'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/4;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 

        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'n-RTD-constM-comp-UDS-s5'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/8;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 

        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'n-RTD-constM-comp-UDS-s6'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/16;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;
        
    case 'n-RTD-constM-comp-UDS-s7'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/32;
        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;
        
    case 'n-RTD-constM-comp-UDS-s8'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/64;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];

        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;
    case 'n-RTD-constM-comp-UDS-s9'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/128;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 

        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        %Vc   = 0*Eg_GaAs+0*[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;
        
    case 'n-RTD-comp'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/2;

        Lc1 = ceil(1*30/a)*a;
        Ls1 = ceil(10/a)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(4/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(10/a)*a;
        Lc2 = ceil(1*30/a)*a; 
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_AlGaAs mc_GaAs mc_AlGaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_a   = [ 0  0 0  0  0 0  0];
        Temp = 300;
        
    case 'RTD2'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/2;

        Lc1 = ceil(30/a)*a;
        Ls1 = ceil(15)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(5/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(15)*a;
        Lc2 = ceil(30/a)*a; 
        
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_AlGaAs mc_GaAs mc_AlGaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_AlGaAs mlh_GaAs mlh_AlGaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_a   = [Ndb 0 0 Ndw 0 0 Ndb];
        N_d   = [ 0  0 0  0  0 0  0];
        Temp = 300;

    case 'RTD-PN'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a/2;

        Lc1 = ceil(30/a)*a;
        Ls1 = ceil(15)*a;
        Lb1 = ceil(3/a)*a;
        Lw  = ceil(5/a)*a;
        Lb2 = ceil(3/a)*a;
        Ls2 = ceil(15)*a;
        Lc2 = ceil(30/a)*a; 
        
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Ls1 Lb1 Lw Lb2 Ls2 Lc2];
        Vc   = Eg_GaAs+[0 0 dEc 0 dEc 0 0];
        Vh  = 0-[0 0 dEv 0 dEv 0 0];
        mc   = [mc_GaAs mc_GaAs mc_AlGaAs mc_GaAs mc_AlGaAs mc_GaAs mc_GaAs];
        mlh  = [mlh_GaAs mlh_GaAs mlh_AlGaAs mlh_GaAs mlh_AlGaAs mlh_GaAs mlh_GaAs];
        mhh = [mhh_GaAs mhh_GaAs mhh_AlGaAs mhh_GaAs mhh_AlGaAs mhh_GaAs mhh_GaAs];
        eps = [eps_GaAs eps_GaAs eps_AlGaAs eps_GaAs eps_AlGaAs eps_GaAs eps_GaAs];
        N_d   = [Ndb 0 0 Ndw 0 0 0];
        N_a   = [ 0  0 0  0 0 0  Ndb];
        Temp = 300;
    case 'LASER-PN'
        
        x_port = 0.3;
        
        a   = 0.5;
        dx  = a;

        Lc1 = ceil(30/a)*a;
        Lw  = ceil(5/a)*a;
        Lc2 = ceil(30/a)*a; 
        
        Ndb = 2E+18;
        Ndw = 0;1E+17;

        Eg_AlGaAs  = 1.424+1.155*x_port+0.37*x_port^2;
        mc_AlGaAs  = 0.063+0.083*x_port;
        mlh_AlGaAs = 0.082+0.068*x_port;
        mhh_AlGaAs = 0.51+0.25*x_port;
        
        eps_AlGaAs = 12.9-2.84*x_port;
        
        Eg_GaAs    = 1.424;
        mc_GaAs    = 0.063;
        mlh_GaAs   = 0.082;
        mhh_GaAs   = 0.51;
        eps_GaAs   = 12.9;
        dEc        = 0.73*(Eg_AlGaAs-Eg_GaAs);
        dEv        = 0.27*(Eg_AlGaAs-Eg_GaAs);
        
        Lx  = [Lc1 Lw Lc2];
        Vc   = Eg_GaAs+[dEc 0 dEc];
        Vh  = 0-[dEv 0 dEv];
        mc   = [mc_AlGaAs mc_GaAs mc_AlGaAs];
        mlh  = [mlh_AlGaAs mlh_GaAs mlh_AlGaAs];
        mhh = [mhh_AlGaAs mhh_GaAs mhh_AlGaAs];
        eps = [eps_AlGaAs eps_GaAs eps_AlGaAs];
        N_d   = [Ndb 0 0];
        N_a   = [0   0 Ndb];
        Temp = 300;
        
    otherwise 
        error('no valid device has been chosen\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAP ON THE 2D COMPUTATIONAL GRID

% Diskretisierung in y-Richtung
if nargin < 3 || isempty(dy)
    dy = dx;
end

% Anzahl der Gitterpunkte in y-Richtung
NyIntervals = max(2, ceil(Ly/dy));
dy = Ly/NyIntervals;
Ny = NyIntervals + 1; % include both endpoints of the periodic cell

% Anzahl der Gitterpunkte pro Schicht in x-Richtung
Nx_layer = round(Lx/dx);

% Gesamtzahl der Gitterpunkte
Nx = sum(Nx_layer) + 1;

% Koordinatenvektoren
x_grid = (0:Nx-1)*dx;
y_grid = (0:Ny-1)*dy;

% 2D-Koordinatenmatrizen
[mat.X, mat.Y] = ndgrid(x_grid, y_grid);
mat.x = x_grid;
mat.y = y_grid;

% Speicher fuer Materialparameter reservieren
mat.Vc  = zeros(Nx, Ny);
mat.Vh  = zeros(Nx, Ny);
mat.mc  = zeros(Nx, Ny);
mat.mlh = zeros(Nx, Ny);
mat.mhh = zeros(Nx, Ny);
mat.eps = zeros(Nx, Ny);

mat.N_d = zeros(Nx, Ny);
mat.N_a = zeros(Nx, Ny);

% Mapping der Materialparameter
IT = 1;

for IL = 1:length(Lx)

    for ix = 1:Nx_layer(IL)

        % Materialparameter werden ueber die gesamte
        % y-Richtung konstant fortgesetzt

        mat.Vc(IT, :)  = Vc(IL);
        mat.Vh(IT, :)  = Vh(IL);

        mat.mc(IT, :)  = mc(IL);
        mat.mlh(IT, :) = mlh(IL);
        mat.mhh(IT, :) = mhh(IL);

        mat.eps(IT, :) = eps(IL);

        mat.N_d(IT, :) = N_d(IL);
        mat.N_a(IT, :) = N_a(IL);

        IT = IT + 1;

    end

end

%% ADDITIONAL DEVICE PARAMETERS

% The last sample belongs to the right lead, not an extra material layer.
fields = {'Vc','Vh','mc','mlh','mhh','eps','N_d','N_a'};
for iField = 1:numel(fields)
    mat.(fields{iField})(end,:) = mat.(fields{iField})(end-1,:);
end
mat.type = device;
mat.V = mat.Vc;
mat.me_x_ch = mc(1);
mat.me_y_ch = mc(1);
mat.me_z_ch = mc(1);
mat.me_x = mat.mc;
mat.me_y = mat.mc;
mat.me_z = mat.mc;
mat.Nd = mat.N_d*1e6; % legacy N_d is cm^-3; full2D uses m^-3
mat.deg_factor = 1;
mat.Vdd = 0;

mat.Temp = Temp;

mat.Lc1 = Lc1;
mat.Lc2 = Lc2;

mat.dx = dx;
mat.dy = dy;

mat.Nx = Nx;
mat.Ny = Ny;

mat.Lx = Lx;
mat.Ly = Ly;

% 1D-Koordinatenvektoren (optional)
mat.x_grid = x_grid;
mat.y_grid = y_grid;
  

