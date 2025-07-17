
function get_parameters(parameters)
    a0=parameters["a0"];
    theta_=parameters["theta_"];
    waa=parameters["waa"];
    wab=parameters["wab"];
    v=parameters["v"];
    valley=parameters["valley"];

    AB_potential=parameters["AB_potential"];
    Phi=parameters["Phi"];
    kappa=parameters["kappa"];#ev, screening of coulomb interaction

    ######################
    theta=theta_/180*pi;
    kd = 4*pi/(3*a0);
    ktheta = 2*kd*sin(theta/2);
    aM = a0/(2*sin(theta/2));
    alpha = wab/(2*v*kd*sin(theta/2));


    return theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa

end




function band_single_k_Helical_trilayer(parameters,k::Vector,nK::Number)
    #diagonalize the continuum model

    theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa=get_parameters(parameters);
    ##########################################

    phi=2*pi/3;

    q1=ktheta*[0, -1];
    q2=ktheta*[sqrt(3)/2, 1/2];
    q3=ktheta*[-sqrt(3)/2, 1/2];
    
    #reciprocal lattice
    b1m=sqrt(3)*ktheta*[0.5, -sqrt(3)/2];
    b2m=sqrt(3)*ktheta*[0.5, sqrt(3)/2];
    
    #momentum shift of three layers
    K1=ktheta*[0, -1];
    K2=ktheta*[0, 0];
    K3=ktheta*[0, 1];
    
    #define tunneling matrices
    n=1;
    T1=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];
    n=2;
    T3=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];
    n=3;
    T2=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];


    T1_12=T1*exp(-im*Phi[1]);
    T2_12=T2*exp(-im*Phi[2]);
    T3_12=T3*exp(-im*Phi[3]);

    T1_23=T1*exp(im*Phi[1]);
    T2_23=T2*exp(im*Phi[2]);
    T3_23=T3*exp(im*Phi[3]);

    if valley==1
        T1_12=conj.(T1_12);
        T2_12=conj.(T2_12);
        T3_12=conj.(T3_12);
    
        T1_23=conj.(T1_23);
        T2_23=conj.(T2_23);
        T3_23=conj.(T3_23);
    end
    
    if valley==-1
        vx_coe=1;
    elseif valley==1
        vx_coe=-1;
    end
    

    ###########################################

    k_set=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set[c1,c2,:]=k+(c1-nK-1)*b1m+(c2-nK-1)*b2m;
        end
    end
    

    layer_set=[1 2 3];#three layers
    sz_set=[1 -1];#A,B,sublattice

    #construct single particle basis
    N_basis=prod([size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_kinetic=zeros((N_basis,N_basis))*im;
    H_U=zeros((N_basis,N_basis))*im;
    H_AB_potential=zeros((N_basis,N_basis))*im;
    # H_kinetic_derivative_x=zeros((N_basis,N_basis))*im;#velocity operator for Chern number
    # H_kinetic_derivative_y=zeros((N_basis,N_basis))*im;#velocity operator for Chern number



    #########################################
    #kinetic energy
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            k=k_set[c1a,c1b,:];
            
            for c2 in eachindex(layer_set)
                for c3 in eachindex(sz_set)
                    if c3==1
                        c31=1;
                        c32=2;
                    elseif c3==2
                        c31=2;
                        c32=1;
                    end
                    
                    ind1=basis[c1a,c1b,c2,c31];
                    ind2=basis[c1a,c1b,c2,c32];


                    if c31==2#sublattice B
                        if c2==1
                            k_=(k-valley*K1);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        elseif c2==2
                            k_=(k-valley*K2);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        elseif c2==3
                            k_=(k-valley*K3);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        end

                    elseif c31==1#sublattice A
                        if c2==1
                            k_=(k-valley*K1);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        elseif c2==2
                            k_=(k-valley*K2);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        elseif c2==3
                            k_=(k-valley*K3);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        end
                    end
                    
                end
            end
        end
    end
    
    #waa and wab in U term
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            for d1a in axes(k_set,1)
                for d1b in axes(k_set,2)
                    for c2 = 1:length(layer_set)-1
                        if c2==1
                            c21=1;
                            c22=2;
                        elseif c2==2
                            c21=2;
                            c22=3;
                        end
                        
                        for c31=1:length(sz_set)#sublattice
                            for c32=1:length(sz_set)#sublattice
                                
                                ind1=basis[c1a,c1b,c21,c31];
                                ind2=basis[d1a,d1b,c22,c32];

                                if c21==1
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K1*valley;
                                elseif c21==2
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K2*valley;
                                elseif c21==3
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K3*valley;
                                end

                                if c22==1
                                    k=k_set[d1a,d1b,:];
                                    k=k-K1*valley;
                                elseif c22==2
                                    k=k_set[d1a,d1b,:];
                                    k=k-K2*valley;
                                elseif c22==3
                                    k=k_set[d1a,d1b,:];
                                    k=k-K3*valley;
                                end

                                if (c21==1)&&(c22==2)
                                    if (sum(abs.(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T1_12[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T2_12[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T3_12[c31,c32];
                                    end
                                elseif (c21==2)&&(c22==3)
                                    if (sum(abs.(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T1_23[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T2_23[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T3_23[c31,c32];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    #sublattice potential to split degenerate bands
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            for c2 in eachindex(layer_set)
                for c3 in eachindex(sz_set) 

                    ind1=basis[c1a,c1b,c2,c3];
                    ind2=basis[c1a,c1b,c2,c3];
                    H_AB_potential[ind1,ind2]=H_AB_potential[ind1,ind2]+AB_potential[c2,c3];
                end
            end
        end
    end
    

    ################################
    H=H_kinetic+H_U+H_U'+H_AB_potential;
    @assert norm(H-H')/norm(H)<1e-12;
    H=(H+H')/2;
    eu,ev=eigen(H);
    return H,eu,ev,b1m,b2m,basis,k_set
end


function Helical_trilayer_kinetic(parameters,k::Vector,nK::Number)


    theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa=get_parameters(parameters);
    ##########################################

    phi=2*pi/3;

    q1=ktheta*[0, -1];
    q2=ktheta*[sqrt(3)/2, 1/2];
    q3=ktheta*[-sqrt(3)/2, 1/2];
    
    b1m=sqrt(3)*ktheta*[0.5, -sqrt(3)/2];
    b2m=sqrt(3)*ktheta*[0.5, sqrt(3)/2];
    
    
    K1=ktheta*[0, -1];
    K2=ktheta*[0, 0];
    K3=ktheta*[0, 1];
    
    if valley==-1
        vx_coe=1;
    elseif valley==1
        vx_coe=-1;
    end

    ###########################################

    k_set=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set[c1,c2,:]=k+(c1-nK-1)*b1m+(c2-nK-1)*b2m;
        end
    end
    

    layer_set=[1 2 3];
    sz_set=[1 -1];#sublattice

    #construct single particle basis
    N_basis=prod([size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_kinetic=zeros((N_basis,N_basis))*im;

    #########################################
    #kinetic energy
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            k=k_set[c1a,c1b,:];
            
            for c2 in eachindex(layer_set)
                for c3 in eachindex(sz_set)
                    if c3==1
                        c31=1;
                        c32=2;
                    elseif c3==2
                        c31=2;
                        c32=1;
                    end
                    
                    ind1=basis[c1a,c1b,c2,c31];
                    ind2=basis[c1a,c1b,c2,c32];


                    if c31==2#sublattice B
                        if c2==1
                            k_=(k-valley*K1);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        elseif c2==2
                            k_=(k-valley*K2);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        elseif c2==3
                            k_=(k-valley*K3);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]-im*k_[2]);
                        end

                    elseif c31==1#sublattice A
                        if c2==1
                            k_=(k-valley*K1);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        elseif c2==2
                            k_=(k-valley*K2);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        elseif c2==3
                            k_=(k-valley*K3);
                            H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(vx_coe*k_[1]+im*k_[2]);
                        end
                    end
                    
                end
            end
        end
    end
    
    

    ################################
    @assert norm(H_kinetic-H_kinetic')/norm(H_kinetic)<1e-12;
    H_kinetic=(H_kinetic+H_kinetic')/2;
    return H_kinetic,b1m,b2m,basis,k_set
end

function Helical_trilayer_potential(parameters,k::Vector,nK::Number)


    theta_,theta,aM,ktheta,waa,wab,v,valley,AB_potential,Phi,kappa=get_parameters(parameters);
    ##########################################

    phi=2*pi/3;

    q1=ktheta*[0, -1];
    q2=ktheta*[sqrt(3)/2, 1/2];
    q3=ktheta*[-sqrt(3)/2, 1/2];
    
    b1m=sqrt(3)*ktheta*[0.5, -sqrt(3)/2];
    b2m=sqrt(3)*ktheta*[0.5, sqrt(3)/2];
    
    
    K1=ktheta*[0, -1];
    K2=ktheta*[0, 0];
    K3=ktheta*[0, 1];
    
    n=1;
    T1=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];
    n=2;
    T3=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];
    n=3;
    T2=[waa wab*exp(-im*phi*(n-1));wab*exp(im*phi*(n-1)) waa];


    T1_12=T1*exp(-im*Phi[1]);
    T2_12=T2*exp(-im*Phi[2]);
    T3_12=T3*exp(-im*Phi[3]);

    T1_23=T1*exp(im*Phi[1]);
    T2_23=T2*exp(im*Phi[2]);
    T3_23=T3*exp(im*Phi[3]);

    if valley==1
        T1_12=conj.(T1_12);
        T2_12=conj.(T2_12);
        T3_12=conj.(T3_12);
    
        T1_23=conj.(T1_23);
        T2_23=conj.(T2_23);
        T3_23=conj.(T3_23);
    end
    
    if valley==-1
        vx_coe=1;
    elseif valley==1
        vx_coe=-1;
    end

    ###########################################

    k_set=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set[c1,c2,:]=k+(c1-nK-1)*b1m+(c2-nK-1)*b2m;
        end
    end
    

    layer_set=[1 2 3];
    sz_set=[1 -1];#sublattice

    #construct single particle basis
    N_basis=prod([size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set,1),size(k_set,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    

    H_U=zeros((N_basis,N_basis))*im;
    H_AB_potential=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_x=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_y=zeros((N_basis,N_basis))*im;



    #########################################
    #kinetic energy
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            k=k_set[c1a,c1b,:];
            
            for c2 in eachindex(layer_set)
                for c3 in eachindex(sz_set)
                    if c3==1
                        c31=1;
                        c32=2;
                    elseif c3==2
                        c31=2;
                        c32=1;
                    end
                    
                    ind1=basis[c1a,c1b,c2,c31];
                    ind2=basis[c1a,c1b,c2,c32];


                    if c31==2#sublattice B
                        if c2==1
                            k_=(k-valley*K1);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]-im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(-im);
                        elseif c2==2
                            k_=(k-valley*K2);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]-im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(-im);
                        elseif c2==3
                            k_=(k-valley*K3);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]-im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(-im);
                        end

                    elseif c31==1#sublattice A
                        if c2==1
                            k_=(k-valley*K1);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]+im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(im);
                        elseif c2==2
                            k_=(k-valley*K2);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]+im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(im);
                        elseif c2==3
                            k_=(k-valley*K3);
                            # H_kinetic[ind1,ind2]=H_kinetic[ind1,ind2]-v*(k_[1]+im*k_[2]);
                            H_kinetic_derivative_x[ind1,ind2]=H_kinetic_derivative_x[ind1,ind2]-vx_coe*v*(1);
                            H_kinetic_derivative_y[ind1,ind2]=H_kinetic_derivative_y[ind1,ind2]-v*(im);
                        end
                    end
                    
                end
            end
        end
    end
    
    #waa and wab in U term
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            for d1a in axes(k_set,1)
                for d1b in axes(k_set,2)
                    for c2 = 1:length(layer_set)-1
                        if c2==1
                            c21=1;
                            c22=2;
                        elseif c2==2
                            c21=2;
                            c22=3;
                        end
                        
                        for c31=1:length(sz_set)#sublattice
                            for c32=1:length(sz_set)#sublattice
                                
                                ind1=basis[c1a,c1b,c21,c31];
                                ind2=basis[d1a,d1b,c22,c32];

                                if c21==1
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K1*valley;
                                elseif c21==2
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K2*valley;
                                elseif c21==3
                                    kp=k_set[c1a,c1b,:];
                                    kp=kp-K3*valley;
                                end

                                if c22==1
                                    k=k_set[d1a,d1b,:];
                                    k=k-K1*valley;
                                elseif c22==2
                                    k=k_set[d1a,d1b,:];
                                    k=k-K2*valley;
                                elseif c22==3
                                    k=k_set[d1a,d1b,:];
                                    k=k-K3*valley;
                                end

                                if (c21==1)&&(c22==2)
                                    if (sum(abs.(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T1_12[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T2_12[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T3_12[c31,c32];
                                    end
                                elseif (c21==2)&&(c22==3)
                                    if (sum(abs.(kp-(k-valley*q1)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T1_23[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q2)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T2_23[c31,c32];
                                    elseif (sum(abs.(kp-(k-valley*q3)))<(1e-10)*norm(q1)) 
                                        H_U[ind1,ind2]=T3_23[c31,c32];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    #sublattice potential
    for c1a in axes(k_set,1)
        for c1b in axes(k_set,2)
            for c2 in eachindex(layer_set)
                for c3 in eachindex(sz_set) 

                    ind1=basis[c1a,c1b,c2,c3];
                    ind2=basis[c1a,c1b,c2,c3];
                    H_AB_potential[ind1,ind2]=H_AB_potential[ind1,ind2]+AB_potential[c2,c3];
                end
            end
        end
    end
    

    ################################
    Hpotential=H_U+H_U'+H_AB_potential;
    @assert norm(Hpotential-Hpotential')/norm(Hpotential)<1e-12;

    return Hpotential,H_kinetic_derivative_x,H_kinetic_derivative_y,b1m,b2m,basis,k_set
end



function get_bandstate(parameters,band_model,N,nK,nb,g1,g3,klabel_set,k1_shift::Vector,k2_shift::Vector,flux_type)

    momentum_set=zeros(size(klabel_set));
    band_state=zeros((N,2*nK+1,2*nK+1,3,2))*im;#k, Kx,Ky,layer,sublattice
    E_band=zeros(N);
    for cc in range(1,size(momentum_set,1))
        momentum_set[cc,:]=klabel_set[cc,1]*g1+klabel_set[cc,2]*g3;
        if flux_type=="spin_independent"
            momentum_shifted=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)-(k2_shift[1]*g1+k2_shift[2]*g3);
            # momentum_shifted_up=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)-(k2_shift[1]*g1+k2_shift[2]*g3);
            # momentum_shifted_dn=momentum_shifted_up;
        elseif flux_type=="spin_dependent"
            # momentum_shifted_up=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)-(k2_shift[1]*g1+k2_shift[2]*g3);
            # momentum_shifted_dn=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)+(k2_shift[1]*g1+k2_shift[2]*g3);
        else
            error("unknown type")
        end
        # momentum_shifted=(momentum_shifted_up,momentum_shifted_dn,);
        #H_0,eu_,ev_,g_set_,kapap_,kapam_,basis_,H_x_,H_y_,g_lattice_=band_model(parameters,momentum_shifted,nK);
        H_0,eu_,ev_,g1,g3,basis_,g_lattice_=band_model(parameters,momentum_shifted,nK);
        

        E_band[cc]=eu_[nb];
        st=reshape(ev_[:,nb],(2*nK+1,2*nK+1,3,2));
        band_state[cc,:,:,:,:]=reshape(st,(size(band_state,2),size(band_state,3),size(band_state,4),size(band_state,5)));
        
    
    end
    return momentum_set,band_state,E_band
end



function get_Fm(N,nK,band_state,momentum_set,klabel_set)
    Fm=zeros((N,N,2*nK+1,2*nK+1))*im;#k1,k2,Kx,Ky
    Fm_label=zeros((N,N,2*nK+1,2*nK+1,2));#k1,k2,Kx,Ky,spin,xy components of momentum
    
    for ind1 in range(1,N)
        for ind2 in range(1,N)
            k1=momentum_set[ind1,:];
            k2=momentum_set[ind2,:];
            q=k1-k2;
            for c1 in range(1,2*nK+1)
                for c2 in range(1,2*nK+1)

                    shift1=c1-(nK+1);
                    shift2=c2-(nK+1);

                    if shift1<=0
                            range_a1=range(1,2*nK+1+shift1);
                            range_a2=range(1-shift1,2*nK+1);
                    elseif shift1>0
                            range_a1=range(1+shift1,2*nK+1);
                            range_a2=range(1,2*nK+1-shift1);
                    end
                    
                    if shift2<=0
                            range_b1=range(1,2*nK+1+shift2);
                            range_b2=range(1-shift2,2*nK+1);
                    elseif shift2>0
                            range_b1=range(1+shift2,2*nK+1);
                            range_b2=range(1,2*nK+1-shift2);
                    end
                    
                    st1=band_state[ind1,range_a1,:,:,:];
                    st1=st1[:,range_b1,:,:];
                    st2=band_state[ind2,range_a2,:,:,:];
                    st2=st2[:,range_b2,:,:];
                    Fm[ind1,ind2,c1,c2]=dot(st1,st2);
                    Fm_label[ind1,ind2,c1,c2,1:2]=klabel_set[ind1,:]-klabel_set[ind2,:]+[shift1,shift2];

                end
            end
        end
    end
    return Fm, Fm_label
end



function search_label(data,klabel,Rank)
    pos=[-1];
    if Rank==3
        (sz1,sz2,sz3)=size(data);
        for (c1,c2) in Iterators.product(1:sz1, 1:sz2)     
            if sum(abs.(data[c1,c2,:]-klabel))<1e-10
                return [c1,c2];
            end
        end
    elseif Rank==4
        (sz1,sz2,sz3,sz4)=size(data);
        for (c1,c2,c3) in Iterators.product(1:sz1, 1:sz2,1:sz3)     
            if sum(abs.(data[c1,c2,c3,:]-klabel))<1e-10
                return [c1,c2,c3];
            end
        end
    end
    return pos    
end



function get_coulomb_potential(parameters,N,nK,aM,g1,g3,klabel_set)
    #coulomb interaction
    q_matrix=zeros((N,2*nK+1,2*nK+1));
    q_matrix_label=zeros((N,2*nK+1,2*nK+1,2));

    k0 = 8.98755e9;
    J_to_meV = 6.24150636e21; #meV
    e_charge = 1.60217663e-19; #coulomb
    epsilon_r = parameters["epsilon"];

    println("epsilon_r="*string(epsilon_r));
    Area = sqrt(3)/2*N*aM^2;
    function coulomb_potential(q,parameters,g1)
        if parameters["screen_type"]=="No"
            return 1/norm(q);
        elseif parameters["screen_type"]=="yukawa"
            return 1/sqrt(q[1]^2+q[2]^2 +(kappa*norm(g1))^2);
            #sqrt(qQ[1]^2+qQ[2]^2+(kappa*norm(g1))^2)
        elseif parameters["screen_type"]=="two_gate"
            return tanh(norm(q)*parameters["d_sc"]*1e-9)*1/norm(q);
        elseif parameters["screen_type"]=="single_gate"
            return (1-exp(-2*norm(q)*parameters["d_sc"]*1e-9))*1/norm(q);
        end
    end
    for c1 in range(1,N)
        for c2 in range(1,2*nK+1)
            for c3 in range(1,2*nK+1)
                q=klabel_set[c1,:];
                Q=[c2-1-nK, c3-1-nK];
                if (sum(abs.(q))==0)&(sum(abs.(Q))==0)
                    q_matrix[c1,c2,c3]=0;
                    q_matrix_label[c1,c2,c3,1:2]=[0,0];
                else
                    qQ=q+Q;
                    q_matrix_label[c1,c2,c3,1:2]=(qQ);
                    qQ=qQ[1]*g1+qQ[2]*g3;
                    #q_matrix[c1,c2,c3] = 2*pi*e_charge^2/(epsilon_r*sqrt(qQ[1]^2+qQ[2]^2+(kappa*norm(g1))^2))*J_to_meV/Area*k0;
                    q_matrix[c1,c2,c3] = 2*pi*e_charge^2/(epsilon_r)*J_to_meV/Area*k0*coulomb_potential(qQ,parameters,g1);
                end
            end
        end
    end
    
    # if parameters["screen_type"]=="two_gate"
    #     @show norm(g1)*parameters["d_sc"]*1e-9
    #     @show tanh(norm(g1)*parameters["d_sc"]*1e-9)
    # end

    V_typical=2*pi*e_charge^2/(epsilon_r*sqrt((norm(g1))^2))*J_to_meV/Area*k0;
    
    return q_matrix,q_matrix_label,V_typical

end




function get_V_k1234_spinup(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234)
    #<k1sigma,k2sigma'|V|k4sigma,k3sigma'>
    Lk=size(k1234,1);
    V_k1234=zeros(Lk)*im;#k1234
    
    for c1 in 1:Lk
        (ind_k1,ind_k2,ind_k3,ind_k4)=k1234[c1,:];
        elem=0;
        for (c2,c3) in Iterators.product(1:2*nK+1, 1:2*nK+1)
            # println(klabel_set[ind_k1,:])
            # println(klabel_set[ind_k4,:])
            dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
            g_label=[c2-1-nK,c3-1-nK];
            q_label=dk_label+g_label;
    
            Fma=Fm_label[ind_k1,ind_k4,:,:,:];
            # posa=search_label(Fma,q_label,3);
            posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
            if length(posa)>0
                posa=posa[1]
            else
                continue
                #posa=[-1];
            end
            Fmb=Fm_label[ind_k2,ind_k3,:,:,:];
            #posb=search_label(Fmb,-q_label,3);
            posb=findall(x->x<1e-10, abs.(Fmb[:,:,1].+q_label[1])+abs.(Fmb[:,:,2].+q_label[2]))
            if length(posb)>0
                posb=posb[1]
            else
                continue
                #posb=[-1];
            end

            #posc=search_label(q_matrix_label,q_label,4);
            posc=findall(x->x<1e-10, abs.(q_matrix_label[:,:,:,1].-q_label[1])+abs.(q_matrix_label[:,:,:,2].-q_label[2]))
            if length(posc)>0
                posc=posc[1]
            else
                continue
                #posc=[-1];
            end


                
            coe1=Fm[ind_k1,ind_k4,posa[1],posa[2]];
            coe2=Fm[ind_k2,ind_k3,posb[1],posb[2]];
            coe3=q_matrix[posc[1],posc[2],posc[3]];
            
            V_k1234[c1]=V_k1234[c1]+coe1*coe2*coe3/2;

        end

        
    end
    return V_k1234
end

 


function get_V_k1234_spinup_single_q(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_label,k1234)

    #cdag_{k1} c_{k1+q}  cdag_{k2}  c_{k2-q}
    Lk=size(k1234,1);
    V_k1234=zeros(Lk)*im;#k1234
    
    for c1 in 1:Lk
        (ind_k1,ind_k4,ind_k2,ind_k3)=k1234[c1,:];

        #for (c2,c3) in Iterators.product(1:2*nK+1, 1:2*nK+1)
            # println(klabel_set[ind_k1,:])
            # println(klabel_set[ind_k4,:])
            dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];

    
            Fma=Fm_label[ind_k1,ind_k4,:,:,:];
            # posa=search_label(Fma,q_label,3);
            # @show Fma[1:2,1:2,:]
            # @show q_label
            posa=findall(x->x<1e-10, abs.(Fma[:,:,1].+q_label[1])+abs.(Fma[:,:,2].+q_label[2]))
            # @show posa
            if length(posa)>0
                posa=posa[1]
            else
                continue
                #posa=[-1];
            end
            Fmb=Fm_label[ind_k2,ind_k3,:,:,:];
            #posb=search_label(Fmb,-q_label,3);
            posb=findall(x->x<1e-10, abs.(Fmb[:,:,1].-q_label[1])+abs.(Fmb[:,:,2].-q_label[2]))
            if length(posb)>0
                posb=posb[1]
            else
                continue
                #posb=[-1];
            end

            #posc=search_label(q_matrix_label,q_label,4);
            posc=findall(x->x<1e-10, abs.(q_matrix_label[:,:,:,1].-q_label[1])+abs.(q_matrix_label[:,:,:,2].-q_label[2]))
            if length(posc)>0
                posc=posc[1]
            else
                continue
                #posc=[-1];
            end


                
            coe1=Fm[ind_k1,ind_k4,posa[1],posa[2]];
            coe2=Fm[ind_k2,ind_k3,posb[1],posb[2]];
            # coe3=q_matrix[posc[1],posc[2],posc[3]];
            
            # V_k1234[c1]=V_k1234[c1]+coe1*coe2*coe3/2;
            V_k1234[c1]=V_k1234[c1]+coe1*coe2;

        #end

        
    end
    return V_k1234
end



