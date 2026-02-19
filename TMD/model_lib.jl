
function band_single_k(theta_::Number,k::Tuple,nK::Number)
    # this is Liang's parameter set
    theta=theta_/180*pi;
    a0=3.52e-10;
    w = 13.3;
    V= 11.2;
    phi=-91/180*pi;
    phi_set=[-phi,phi];

    me = 9.10938356e-31; # kg
    m = 0.62*me;
    hbar = 1.054571817e-34;
    J_to_meV = 6.24150636e21; #meV
    mass_coe = hbar^2/(2*m)*J_to_meV;


    aM=a0/(2*sin(theta/2));
    #aM=a0/(theta);#approximation
    g_set=zeros((6,2));
    # plot(0,0,'ro');hold on;
    for cc in range(1,size(g_set,1))
        vv=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)]
        g_set[cc,:]=4*pi/(sqrt(3)*aM)*vv;
    end
    g1=g_set[1,:];
    g2=g_set[2,:];
    g3=g_set[3,:];
    g4=g_set[4,:];
    g5=g_set[5,:];
    g6=g_set[6,:];
    kapa2=(g1+g6)/3;#kapa-
    kapa1=(g1+g2)/3;#kapa+

    k_set_up=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_up[c1,c2,:]=k[1]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end
    k_set_dn=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_dn[c1,c2,:]=k[2]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end
    

    layer_set=[1,2];
    sz_set=[1,-1];

    #construct single particle basis
    N_basis=prod([size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_kinetic=zeros((N_basis,N_basis))*im;
    H_V=zeros((N_basis,N_basis))*im;
    H_t=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_x=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_y=zeros((N_basis,N_basis))*im;


    #kinetic energy
    for c1a in range(1,size(k_set_up,1))
        for c1b in range(1,size(k_set_up,2))

            for c2 in range(1,length(layer_set))
                for c3 in range(1,length(sz_set))
                    if c3==1#spin up
                        k=k_set_up[c1a,c1b,:];
                    elseif c3==2#spin dn
                        k=k_set_dn[c1a,c1b,:];
                    end
                    k=k[:]


                    ind=Int64.(basis[c1a,c1b,c2,c3]);
                    # println(k)
                    # println(ind)
                    if c2==1
                        kapa=kapa1;
                    elseif c2==2
                        kapa=kapa2;
                    end
                    

                    if c3==1#spin up
                        # println(sum((k-kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k-kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]-kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]-kapa[2]))*mass_coe;
                    elseif c3==2#spin down
                        # println(sum((k+kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k+kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]+kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]+kapa[2]))*mass_coe;
                    end
                end
            end
        end
    end
                    


    #V term
    for c3 in range(1,length(sz_set))
        if c3==1
            k_set=k_set_up;
        elseif c3==2
            k_set=k_set_dn;
        end

        for c1a in range(1,size(k_set,1))
            for c1b in range(1,size(k_set,2))
                kp=k_set[c1a,c1b,:];
                kp=kp[:];
                for d1a in range(1,size(k_set,1))
                    for d1b in range(1,size(k_set,2))
                        k=k_set[d1a,d1b,:];
                        k=k[:];
                        if sum(abs.(kp-(k+g1)))<(1e-10)*norm(g1)
                            jj=1;
                        elseif sum(abs.(kp-(k+g2)))<(1e-10)*norm(g1)
                            jj=2;
                        elseif sum(abs.(kp-(k+g3)))<(1e-10)*norm(g1)
                            jj=3;
                        elseif sum(abs.(kp-(k+g4)))<(1e-10)*norm(g1)
                            jj=4;
                        elseif sum(abs.(kp-(k+g5)))<(1e-10)*norm(g1)
                            jj=5;
                        elseif sum(abs.(kp-(k+g6)))<(1e-10)*norm(g1)
                            jj=6;
                        else
                            continue
                        end
                        
                        for c2 in range(1,length(layer_set))
                            
                                ind1=Int64.(basis[c1a,c1b,c2,c3]);
                                ind2=Int64.(basis[d1a,d1b,c2,c3]);
                                H_V[ind1,ind2]=-V*exp(-(im)*((-1)^(jj+1))*phi_set[c2]);#add a minus sign to get the sam result as Yang Peng
                            
                        end
                    end
                end
            end
        end
    end



    #t term for spin up
    k_set=k_set_up;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];

                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==1)&(d2==2)
                                if sum(abs.(kp-(k)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+g2)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+g3)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    coe=0;
                                end
                                
                            elseif (c2==2)&(d2==1)
                                if sum(abs.(kp-(k)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-g2)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-g3)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    coe=0;
                                end
                            else
                                continue
                            end
    
                            c3=1;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=w*coe;
                        end
                    end
                end
            end
        end
    end



    #t term for spin down
    k_set=k_set_dn;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];
                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==2)&(d2==1)
                                if sum(abs.(kp-(k)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+g2)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+g3)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    coe=0;
                                end
                                
                            elseif (c2==1)&(d2==2)
                                if sum(abs.(kp-(k)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-g2)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-g3)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    coe=0;
                                end
                            else
                                continue
                            end
                            
    
                            c3=2;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=w*coe;
                        end
                    end
                end
            end
        end
    end

    H=H_kinetic+H_V+H_t;
    eu,ev=eigen(H);

    return H,eu,ev,g_set,kapa1,kapa2,basis,H_kinetic_derivative_x,H_kinetic_derivative_y,(k_set_up,k_set_dn,)
end


function band_single_k_2nd_harmonic(theta_::Number,k::Tuple,nK::Number)

    theta=theta_/180*pi;
    a0=3.52e-10;


    V=-17.5;
    V2=10.5;
    w=-6.5;
    w2=11.75;
    phi=-58.49/180*pi;
    
    # V=2.4;
    # V2=1.0;
    # w=-5.8;
    # w2=2.8;
    # phi=-90/180*pi;


    phi_set=[-phi,phi];

    me = 9.10938356e-31; # kg
    m = 0.62*me;
    hbar = 1.054571817e-34;
    J_to_meV = 6.24150636e21; #meV
    mass_coe = hbar^2/(2*m)*J_to_meV;


    aM=a0/(2*sin(theta/2));
    #aM=a0/(theta);#approximation
    g_set=zeros((6,2));
    # plot(0,0,'ro');hold on;
    for cc in range(1,size(g_set,1))
        vv=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)]
        g_set[cc,:]=4*pi/(sqrt(3)*aM)*vv;
    end
    g1=g_set[1,:];
    g2=g_set[2,:];
    g3=g_set[3,:];
    g4=g_set[4,:];
    g5=g_set[5,:];
    g6=g_set[6,:];
    kapa2=(g1+g6)/3;#kapa-
    kapa1=(g1+g2)/3;#kapa+


    C3_rotation=[cos(2*pi/3) -sin(2*pi/3);sin(2*pi/3) cos(2*pi/3)];
    q1=(g2-g6)/3;
    q2=(C3_rotation*q1);
    q3=(C3_rotation*C3_rotation*q1);
    
    g21=(g1+g2);
    C3_g21=(C3_rotation*g21);
    C3C3_g21=(C3_rotation*C3_rotation*g21);
    
    q21=(q1+g1);
    C3_q21=(C3_rotation*q21);
    C3C3_q21=(C3_rotation*C3_rotation*q21);




    k_set_up=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_up[c1,c2,:]=k[1]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end
    k_set_dn=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_dn[c1,c2,:]=k[2]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end
    

    layer_set=[1,2];
    sz_set=[1,-1];

    #construct single particle basis
    N_basis=prod([size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_kinetic=zeros((N_basis,N_basis))*im;
    H_V=zeros((N_basis,N_basis))*im;
    H_t=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_x=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_y=zeros((N_basis,N_basis))*im;


    #kinetic energy
    for c1a in range(1,size(k_set_up,1))
        for c1b in range(1,size(k_set_up,2))

            for c2 in range(1,length(layer_set))
                for c3 in range(1,length(sz_set))
                    if c3==1#spin up
                        k=k_set_up[c1a,c1b,:];
                    elseif c3==2#spin dn
                        k=k_set_dn[c1a,c1b,:];
                    end
                    k=k[:]


                    ind=Int64.(basis[c1a,c1b,c2,c3]);
                    # println(k)
                    # println(ind)
                    if c2==1
                        kapa=kapa1;
                    elseif c2==2
                        kapa=kapa2;
                    end
                    

                    if c3==1#spin up
                        # println(sum((k-kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k-kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]-kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]-kapa[2]))*mass_coe;
                    elseif c3==2#spin down
                        # println(sum((k+kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k+kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]+kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]+kapa[2]))*mass_coe;
                    end
                end
            end
        end
    end
                    


    #V term: first harmonic
    for c3 in range(1,length(sz_set))
        if c3==1
            k_set=k_set_up;
        elseif c3==2
            k_set=k_set_dn;
        end

        for c1a in range(1,size(k_set,1))
            for c1b in range(1,size(k_set,2))
                kp=k_set[c1a,c1b,:];
                kp=kp[:];
                for d1a in range(1,size(k_set,1))
                    for d1b in range(1,size(k_set,2))
                        k=k_set[d1a,d1b,:];
                        k=k[:];
                        if sum(abs.(kp-(k+g1)))<(1e-10)*norm(g1)
                            jj=1;
                        elseif sum(abs.(kp-(k+g2)))<(1e-10)*norm(g1)
                            jj=2;
                        elseif sum(abs.(kp-(k+g3)))<(1e-10)*norm(g1)
                            jj=3;
                        elseif sum(abs.(kp-(k+g4)))<(1e-10)*norm(g1)
                            jj=4;
                        elseif sum(abs.(kp-(k+g5)))<(1e-10)*norm(g1)
                            jj=5;
                        elseif sum(abs.(kp-(k+g6)))<(1e-10)*norm(g1)
                            jj=6;
                        else
                            continue
                        end
                        
                        for c2 in range(1,length(layer_set))
                            
                                ind1=Int64.(basis[c1a,c1b,c2,c3]);
                                ind2=Int64.(basis[d1a,d1b,c2,c3]);
                                H_V[ind1,ind2]=H_V[ind1,ind2]+V*exp(-(im)*((-1)^(jj+1))*phi_set[c2]);#add a minus sign to get the sam result as Yang Peng
                            
                        end
                    end
                end
            end
        end
    end

    #V2 term: second harmonic
    for c3 in range(1,length(sz_set))
        if c3==1
            k_set=k_set_up;
        elseif c3==2
            k_set=k_set_dn;
        end

        for c1a in range(1,size(k_set,1))
            for c1b in range(1,size(k_set,2))
                kp=k_set[c1a,c1b,:];
                kp=kp[:];
                for d1a in range(1,size(k_set,1))
                    for d1b in range(1,size(k_set,2))
                        k=k_set[d1a,d1b,:];
                        k=k[:];
                        if sum(abs.(kp-(k-g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k-C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k-C3C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+C3C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        else
                            continue
                        end
                        
                        for c2 in range(1,length(layer_set))
                            
                                ind1=Int64.(basis[c1a,c1b,c2,c3]);
                                ind2=Int64.(basis[d1a,d1b,c2,c3]);
                                H_V[ind1,ind2]=H_V[ind1,ind2]+V2*coe;#add a minus sign to get the sam result as Yang Peng
                            
                        end
                    end
                end
            end
        end
    end



    #t term for spin up: first harmonic
    k_set=k_set_up;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];

                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==1)&(d2==2)
                                if sum(abs.(kp-(k-q1+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q2+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q3+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==2)&(d2==1)
                                if sum(abs.(kp-(k+q1-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q2-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q3-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
    
                            c3=1;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w*coe;
                        end
                    end
                end
            end
        end
    end

    #t term for spin down: first harmonic
    k_set=k_set_dn;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];
                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==2)&(d2==1)
                                if sum(abs.(kp-(k-q1+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q2+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q3+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==1)&(d2==2)
                                if sum(abs.(kp-(k+q1-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q2-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q3-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
                            
    
                            c3=2;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w*coe;
                        end
                    end
                end
            end
        end
    end


    #t term for spin up: second harmonic
    k_set=k_set_up;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];

                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==1)&(d2==2)
                                if sum(abs.(kp-(k-q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==2)&(d2==1)
                                if sum(abs.(kp-(k+q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
    
                            c3=1;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w2*coe;
                        end
                    end
                end
            end
        end
    end

    #t term for spin down: second harmonic
    k_set=k_set_dn;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];
                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==2)&(d2==1)
                                if sum(abs.(kp-(k-q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==1)&(d2==2)
                                if sum(abs.(kp-(k+q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
                            
    
                            c3=2;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w2*coe;
                        end
                    end
                end
            end
        end
    end

    H=H_kinetic+H_V+H_t;
    eu,ev=eigen(H);

    return H,eu,ev,g_set,kapa1,kapa2,basis,H_kinetic_derivative_x,H_kinetic_derivative_y,(k_set_up,k_set_dn,)
end

function band_single_k_2nd_harmonic_kinetic(k,nK::Number,parameters)
    V=parameters["V"];
    V2=parameters["V2"];
    w=parameters["w"];
    w2=parameters["w2"];
    phi=parameters["phi"];
    theta_=parameters["theta_"];

    theta=theta_/180*pi;
    a0=3.52e-10;



    phi_set=[-phi,phi];

    me = 9.10938356e-31; # kg
    m = 0.62*me;
    hbar = 1.054571817e-34;
    J_to_meV = 6.24150636e21; #meV
    mass_coe = hbar^2/(2*m)*J_to_meV;


    aM=a0/(2*sin(theta/2));
    #aM=a0/(theta);#approximation
    g_set=zeros((6,2));
    # plot(0,0,'ro');hold on;
    for cc in range(1,size(g_set,1))
        vv=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)]
        g_set[cc,:]=4*pi/(sqrt(3)*aM)*vv;
    end
    g1=g_set[1,:];
    g2=g_set[2,:];
    g3=g_set[3,:];
    g4=g_set[4,:];
    g5=g_set[5,:];
    g6=g_set[6,:];
    kapa2=(g1+g6)/3;#kapa-
    kapa1=(g1+g2)/3;#kapa+


    C3_rotation=[cos(2*pi/3) -sin(2*pi/3);sin(2*pi/3) cos(2*pi/3)];
    q1=(g2-g6)/3;
    q2=(C3_rotation*q1);
    q3=(C3_rotation*C3_rotation*q1);
    
    g21=(g1+g2);
    C3_g21=(C3_rotation*g21);
    C3C3_g21=(C3_rotation*C3_rotation*g21);
    
    q21=(q1+g1);
    C3_q21=(C3_rotation*q21);
    C3C3_q21=(C3_rotation*C3_rotation*q21);




    k_set_up=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_up[c1,c2,:]=k[1]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end

    

    layer_set=[1,2];
    sz_set=[1];

    #construct single particle basis
    N_basis=prod([size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_kinetic=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_x=zeros((N_basis,N_basis))*im;
    H_kinetic_derivative_y=zeros((N_basis,N_basis))*im;


    #kinetic energy
    for c1a in range(1,size(k_set_up,1))
        for c1b in range(1,size(k_set_up,2))

            for c2 in range(1,length(layer_set))
                for c3 in range(1,length(sz_set))
                    if c3==1#spin up
                        k=k_set_up[c1a,c1b,:];
                    elseif c3==2#spin dn
                        k=k_set_dn[c1a,c1b,:];
                    end
                    k=k[:]


                    ind=Int64.(basis[c1a,c1b,c2,c3]);
                    # println(k)
                    # println(ind)
                    if c2==1
                        kapa=kapa1;
                    elseif c2==2
                        kapa=kapa2;
                    end
                    

                    if c3==1#spin up
                        # println(sum((k-kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k-kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]-kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]-kapa[2]))*mass_coe;
                    elseif c3==2#spin down
                        # println(sum((k+kapa)^2)*mass_coe)
                        H_kinetic[ind,ind]=sum((k+kapa).^2)*mass_coe;
                        H_kinetic_derivative_x[ind,ind]=2*((k[1]+kapa[1]))*mass_coe;
                        H_kinetic_derivative_y[ind,ind]=2*((k[2]+kapa[2]))*mass_coe;
                    end
                end
            end
        end
    end
                    

    return H_kinetic,H_kinetic_derivative_x,H_kinetic_derivative_y
end

function band_single_k_2nd_harmonic_potential(k::Tuple,nK::Number,parameters)
    V=parameters["V"];
    V2=parameters["V2"];
    w=parameters["w"];
    w2=parameters["w2"];
    phi=parameters["phi"];
    theta_=parameters["theta_"];

    theta=theta_/180*pi;
    a0=3.52e-10;



    phi_set=[-phi,phi];

    me = 9.10938356e-31; # kg
    m = 0.62*me;
    hbar = 1.054571817e-34;
    J_to_meV = 6.24150636e21; #meV
    mass_coe = hbar^2/(2*m)*J_to_meV;


    aM=a0/(2*sin(theta/2));
    #aM=a0/(theta);#approximation
    g_set=zeros((6,2));
    # plot(0,0,'ro');hold on;
    for cc in range(1,size(g_set,1))
        vv=[cos(pi*(cc-1)/3),sin(pi*(cc-1)/3)]
        g_set[cc,:]=4*pi/(sqrt(3)*aM)*vv;
    end
    g1=g_set[1,:];
    g2=g_set[2,:];
    g3=g_set[3,:];
    g4=g_set[4,:];
    g5=g_set[5,:];
    g6=g_set[6,:];
    kapa2=(g1+g6)/3;#kapa-
    kapa1=(g1+g2)/3;#kapa+


    C3_rotation=[cos(2*pi/3) -sin(2*pi/3);sin(2*pi/3) cos(2*pi/3)];
    q1=(g2-g6)/3;
    q2=(C3_rotation*q1);
    q3=(C3_rotation*C3_rotation*q1);
    
    g21=(g1+g2);
    C3_g21=(C3_rotation*g21);
    C3C3_g21=(C3_rotation*C3_rotation*g21);
    
    q21=(q1+g1);
    C3_q21=(C3_rotation*q21);
    C3C3_q21=(C3_rotation*C3_rotation*q21);


    k_set_up=zeros((2*nK+1,2*nK+1,2));
    for c1 in range(1,2*nK+1)
        for c2 in range(1,2*nK+1)
            k_set_up[c1,c2,:]=k[1]+(c1-nK-1)*g1+(c2-nK-1)*g3;
        end
    end

    layer_set=[1,2];
    sz_set=[1];

    #construct single particle basis
    N_basis=prod([size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)]);
    basis=range(1,N_basis);
    basis=reshape(basis,(size(k_set_up,1),size(k_set_up,2),length(layer_set),length(sz_set)));
    basis=Int64.(basis)
    
    H_V=zeros((N_basis,N_basis))*im;
    H_t=zeros((N_basis,N_basis))*im;


    #V term: first harmonic
    for c3 in range(1,length(sz_set))
        if c3==1
            k_set=k_set_up;
        elseif c3==2
            k_set=k_set_dn;
        end

        for c1a in range(1,size(k_set,1))
            for c1b in range(1,size(k_set,2))
                kp=k_set[c1a,c1b,:];
                kp=kp[:];
                for d1a in range(1,size(k_set,1))
                    for d1b in range(1,size(k_set,2))
                        k=k_set[d1a,d1b,:];
                        k=k[:];
                        if sum(abs.(kp-(k+g1)))<(1e-10)*norm(g1)
                            jj=1;
                        elseif sum(abs.(kp-(k+g2)))<(1e-10)*norm(g1)
                            jj=2;
                        elseif sum(abs.(kp-(k+g3)))<(1e-10)*norm(g1)
                            jj=3;
                        elseif sum(abs.(kp-(k+g4)))<(1e-10)*norm(g1)
                            jj=4;
                        elseif sum(abs.(kp-(k+g5)))<(1e-10)*norm(g1)
                            jj=5;
                        elseif sum(abs.(kp-(k+g6)))<(1e-10)*norm(g1)
                            jj=6;
                        else
                            continue
                        end
                        
                        for c2 in range(1,length(layer_set))
                            
                                ind1=Int64.(basis[c1a,c1b,c2,c3]);
                                ind2=Int64.(basis[d1a,d1b,c2,c3]);
                                H_V[ind1,ind2]=H_V[ind1,ind2]+V*exp(-(im)*((-1)^(jj+1))*phi_set[c2]);#add a minus sign to get the sam result as Yang Peng
                            
                        end
                    end
                end
            end
        end
    end

    #V2 term: second harmonic
    for c3 in range(1,length(sz_set))
        if c3==1
            k_set=k_set_up;
        elseif c3==2
            k_set=k_set_dn;
        end

        for c1a in range(1,size(k_set,1))
            for c1b in range(1,size(k_set,2))
                kp=k_set[c1a,c1b,:];
                kp=kp[:];
                for d1a in range(1,size(k_set,1))
                    for d1b in range(1,size(k_set,2))
                        k=k_set[d1a,d1b,:];
                        k=k[:];
                        if sum(abs.(kp-(k-g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k-C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k-C3C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        elseif sum(abs.(kp-(k+C3C3_g21)))<(1e-10)*norm(g1)
                            coe=1;
                        else
                            continue
                        end
                        
                        for c2 in range(1,length(layer_set))
                            
                                ind1=Int64.(basis[c1a,c1b,c2,c3]);
                                ind2=Int64.(basis[d1a,d1b,c2,c3]);
                                H_V[ind1,ind2]=H_V[ind1,ind2]+V2*coe;#add a minus sign to get the sam result as Yang Peng
                            
                        end
                    end
                end
            end
        end
    end



    #t term for spin up: first harmonic
    k_set=k_set_up;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];

                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==1)&(d2==2)
                                if sum(abs.(kp-(k-q1+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q2+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-q3+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==2)&(d2==1)
                                if sum(abs.(kp-(k+q1-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q2-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+q3-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
    
                            c3=1;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w*coe;
                        end
                    end
                end
            end
        end
    end



    #t term for spin up: second harmonic
    k_set=k_set_up;
    for c1a in range(1,size(k_set,1))
        for c1b in range(1,size(k_set,2))
            kp=k_set[c1a,c1b,:];
            kp=kp[:];
            for d1a in range(1,size(k_set,1))
                for d1b in range(1,size(k_set,2))
                    k=k_set[d1a,d1b,:];
                    k=k[:];

                    for c2 in range(1,length(layer_set))
                        for d2 in range(1,length(layer_set))
                            if (c2==1)&(d2==2)
                                if sum(abs.(kp-(k-q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k-C3C3_q21+q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                                
                            elseif (c2==2)&(d2==1)
                                if sum(abs.(kp-(k+q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                elseif sum(abs.(kp-(k+C3C3_q21-q1)))<(1e-10)*norm(g1)
                                    coe=1;
                                else
                                    continue
                                end
                            else
                                continue
                            end
    
                            c3=1;
                            ind1=Int64.(basis[c1a,c1b,c2,c3]);
                            ind2=Int64.(basis[d1a,d1b,d2,c3]);

                            H_t[ind1,ind2]=H_t[ind1,ind2]+w2*coe;
                        end
                    end
                end
            end
        end
    end

    H_potential=H_V+H_t;
    return H_potential
end

function get_bandstate(band_model,N,nK,N_bands,g1,g3,theta_,klabel_set,k1_shift::Vector,k2_shift::Vector,flux_type)

    momentum_set=zeros(size(klabel_set));
    band_state=zeros((N,2*nK+1,2*nK+1,2,2,N_bands))*im;#k, Kx,Ky,layer,spin,band
    E_band=zeros((N,2,N_bands));#Kx,Ky,layer,spin,band
    for cc in range(1,size(momentum_set,1))
        momentum_set[cc,:]=klabel_set[cc,1]*g1+klabel_set[cc,2]*g3;
        if flux_type=="spin_independent"
            momentum_shifted_up=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)-(k2_shift[1]*g1+k2_shift[2]*g3);
            momentum_shifted_dn=momentum_shifted_up;
        elseif flux_type=="spin_dependent"
            momentum_shifted_up=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)-(k2_shift[1]*g1+k2_shift[2]*g3);
            momentum_shifted_dn=momentum_set[cc,:]-(k1_shift[1]*g1+k1_shift[2]*g3)+(k2_shift[1]*g1+k2_shift[2]*g3);
        else
            error("unknown type")
        end
        momentum_shifted=(momentum_shifted_up,momentum_shifted_dn,);
        H_0,eu_,ev_,g_set_,kapap_,kapam_,basis_,H_x_,H_y_,g_lattice_=band_model(theta_,momentum_shifted,nK);
        
        for spin_ in range(1,2)
            ind_=basis_[:,:,:,spin_];
            ind_=ind_[:]
            H_=H_0[ind_,:];
            H_=H_[:,ind_];
            H_=(H_+H_')/2;
            eu_,ev_=eigen(H_);
            
            for nb in range(1,N_bands)
                E_band[cc,spin_,nb]=eu_[nb];
                st=reshape(ev_[:,nb],(2*nK+1,2*nK+1,2));
                if abs(st[nK+1,nK+1,1])>1e-12
                    st=st/((st[nK+1,nK+1,1])/abs(st[nK+1,nK+1,1]));#choose gauge for state vector
                end
                band_state[cc,:,:,:,spin_,nb]=reshape(st,(size(band_state,2),size(band_state,3),size(band_state,4)));
            end
        end
    end
    return momentum_set,band_state,E_band
end







# function get_Fm(N,nK,N_bands,band_state,momentum_set,klabel_set)
#     Fm=zeros((N,N,2*nK+1,2*nK+1,2,N_bands))*im;#k1,k2,Kx,Ky,spin,band
#     Fm_label=zeros((N,N,2*nK+1,2*nK+1,2,2));#k1,k2,Kx,Ky,spin,xy direction
    
#     for ind1 in range(1,N)
#         for ind2 in range(1,N)
#             k1=momentum_set[ind1,:];
#             k2=momentum_set[ind2,:];
#             q=k1-k2;
#             for c1 in range(1,2*nK+1)
#                 for c2 in range(1,2*nK+1)
#                     for spin_ in range(1,2)
#                         for nb in range(1,N_bands)
#                             shift1=c1-(nK+1);
#                             shift2=c2-(nK+1);
    
#                             if shift1<=0
#                                     range_a1=range(1,2*nK+1+shift1);
#                                     range_a2=range(1-shift1,2*nK+1);
#                             elseif shift1>0
#                                     range_a1=range(1+shift1,2*nK+1);
#                                     range_a2=range(1,2*nK+1-shift1);
#                             end
                            
#                             if shift2<=0
#                                     range_b1=range(1,2*nK+1+shift2);
#                                     range_b2=range(1-shift2,2*nK+1);
#                             elseif shift2>0
#                                     range_b1=range(1+shift2,2*nK+1);
#                                     range_b2=range(1,2*nK+1-shift2);
#                             end
                            
#                             st1=band_state[ind1,range_a1,:,:,spin_,nb];
#                             st1=st1[:,range_b1,:];
#                             st2=band_state[ind2,range_a2,:,:,spin_,nb];
#                             st2=st2[:,range_b2,:];
#                             Fm[ind1,ind2,c1,c2,spin_,nb]=dot(st1,st2);
#                             Fm_label[ind1,ind2,c1,c2,spin_,1:2]=klabel_set[ind1,:]-klabel_set[ind2,:]+[shift1,shift2];
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return Fm, Fm_label
# end


function get_Fm(N,nK,N_bands,band_state,momentum_set,klabel_set)
    Fm=zeros((N,N,2*nK+1,2*nK+1,2,N_bands,N_bands))*im;#k1,k2,Kx,Ky,spin,band,band
    Fm_label=zeros((N,N,2*nK+1,2*nK+1,2,2));#k1,k2,Kx,Ky,spin,xy direction
    
    for ind1 in range(1,N)
        for ind2 in range(1,N)
            k1=momentum_set[ind1,:];
            k2=momentum_set[ind2,:];
            q=k1-k2;
            for c1 in range(1,2*nK+1)
                for c2 in range(1,2*nK+1)
                    for spin_ in range(1,2)
                        for nb1 in range(1,N_bands)
                            for nb2 in range(1,N_bands)
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
                                
                                st1=band_state[ind1,range_a1,:,:,spin_,nb1];
                                st1=st1[:,range_b1,:];
                                st2=band_state[ind2,range_a2,:,:,spin_,nb2];
                                st2=st2[:,range_b2,:];
                                Fm[ind1,ind2,c1,c2,spin_,nb1,nb2]=dot(st1,st2);
                                Fm_label[ind1,ind2,c1,c2,spin_,1:2]=klabel_set[ind1,:]-klabel_set[ind2,:]+[shift1,shift2];
                            end
                        end
                    end
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



# function search_label(data,klabel,Rank)

#     if Rank==3
#         @tensor km[:]:=ones((size(data,1),size(data,2)))[-1,-2]*klabel[-3]
#         diff=sum(abs.(data-km),dims=Rank);
#         #pos=np.where(diff < 1e-10)
#         pos=findall(x->x<1e-10, diff)        
        
#         if length(pos)>0
#             pos=pos[1];
#             pos=[pos[1],pos[2]];
#         else
#             pos=[-1]
#         end

#     elseif Rank==4
#         @tensor km[:]:=ones((size(data,1),size(data,2),size(data,3)))[-1,-2,-3]*klabel[-4]
#         diff=sum(abs.(data-km),dims=Rank);
#         # pos=np.where(diff < 1e-10)
#         pos=findall(x->x<1e-10, diff)
        
#         if length(pos)>0
#             pos=pos[1];
#             pos=[pos[1],pos[2],pos[3]];
#         else
#             pos=[-1]
#         end
#     end
#     #println(pos)
#     return pos
# end
    




# function get_V_k1234_spinful(nK,nb,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234)
#     #<k1sigma,k2sigma'|V|k4sigma,k3sigma'>
#     V_k1234=zeros((size(k1234,1),2,2))*im;#k1234,sigma,sigma'
    

#     for spin1 in range(1,2)
#         for spin2 in range(1,2)
#             for c1 in range(1,size(k1234,1))#range(10):
#                 ind_k1=k1234[c1,0];
#                 ind_k2=k1234[c1,1];
#                 ind_k3=k1234[c1,2];
#                 ind_k4=k1234[c1,3];
    
#                 elem=0;
#                 for c2 in range(1,2*nK+1)
#                     for c3 in range(1,2*nK+1)
#                         # println(klabel_set[ind_k1,:])
#                         # println(klabel_set[ind_k4,:])
#                         dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
#                         g_label=[c2-1-nK,c3-1-nK];
#                         q_label=dk_label+g_label;
                
#                         Fma=reshape(Fm_label[ind_k1,ind_k4,:,:,spin1,nb,:],(2*nK+1,2*nK+1,2));
#                         posa=search_label(Fma,q_label,3);
#                         Fmb=reshape(Fm_label[ind_k2,ind_k3,:,:,spin2,nb,:],(2*nK+1,2*nK+1,2));
#                         posb=search_label(Fmb,-q_label,3);
    
#                         posc=search_label(q_matrix_label,q_label,4);
                        
#                         if (len(posa)==2)&(len(posb)==2)&(len(posc)==3)
#                             #[posa;posb];
#                             #posc;
                            
#                             coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb];
#                             coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb];
#                             coe3=q_matrix[posc[1],posc[2],posc[3]];
                            
#                             elem=elem+coe1*coe2*coe3;
#                         else
#                             if sum(abs(q_label))<=nK
#                                 q_label
#                                 error("not find")
#                             end
#                         end
#                     end
#                 end
#                 V_k1234[c1,spin1,spin2]=elem/2;    
#             end
#         end
#     end
#     return V_k1234
# end


# function get_V_k1234_spinup(nK,nb,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234)
#     #<k1sigma,k2sigma'|V|k4sigma,k3sigma'>
#     Lk=size(k1234,1);
#     V_k1234=zeros((Lk,2,2))*im;#k1234,sigma,sigma'
    

#     for spin1 in range(1,2)
#         for spin2 in range(1,2)
#             for c1 in range(1,Lk)
#                 (ind_k1,ind_k2,ind_k3,ind_k4)=k1234[c1,:];

    
#                 elem=0;
#                 for (c2,c3) in Iterators.product(1:2*nK+1, 1:2*nK+1)
#                 # for c2 in range(1,2*nK+1)
#                 #     for c3 in range(1,2*nK+1)
#                         # println(klabel_set[ind_k1,:])
#                         # println(klabel_set[ind_k4,:])
#                         dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
#                         g_label=[c2-1-nK,c3-1-nK];
#                         q_label=dk_label+g_label;
                
#                         Fma=reshape(Fm_label[ind_k1,ind_k4,:,:,spin1,nb,:],(2*nK+1,2*nK+1,2));
#                         posa=search_label(Fma,q_label,3);
#                         Fmb=reshape(Fm_label[ind_k2,ind_k3,:,:,spin2,nb,:],(2*nK+1,2*nK+1,2));
#                         posb=search_label(Fmb,-q_label,3);
    
#                         posc=search_label(q_matrix_label,q_label,4);
                        
#                         if (length(posa)==2)&&(length(posb)==2)&&(length(posc)==3)
#                             #[posa;posb];
#                             #posc;
                            
#                             coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb];
#                             coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb];
#                             coe3=q_matrix[posc[1],posc[2],posc[3]];
                            
#                             elem=elem+coe1*coe2*coe3;
#                         else
#                             if sum(abs.(q_label))<=nK
#                                 q_label
#                                 error("not find")
#                             end
#                         end
#                 #     end
#                 # end
#                 end
    
#                 V_k1234[c1,spin1,spin2]=elem/2;
#             end
#         end
#     end
#     return V_k1234
# end

function get_V_k1234_spinup(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234)
    #<k1sigma,k2sigma'|V|k4sigma,k3sigma'>
    Lk=size(k1234,1);
    V_k1234=zeros(Lk,2,2,2)*im;#k1234,sigma,sigma'
    
    for (spin1, spin2, c1) in Iterators.product(1:1, 1:1, 1:Lk)
        (ind_k1,ind_k2,ind_k3,ind_k4)=k1234[c1,:];
        elem=0;
        for (c2,c3) in Iterators.product(1:2*nK+1, 1:2*nK+1)
            # println(klabel_set[ind_k1,:])
            # println(klabel_set[ind_k4,:])
            dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
            g_label=[c2-1-nK,c3-1-nK];
            q_label=dk_label+g_label;
    
            Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
            # posa=search_label(Fma,q_label,3);
            posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
            if length(posa)>0
                posa=posa[1]
            else
                continue
                #posa=[-1];
            end
            Fmb=Fm_label[ind_k2,ind_k3,:,:,spin2,:];
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

            for nb1 in range(1,2)
                
                coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb1];
                coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb1,nb1];
                coe3=q_matrix[posc[1],posc[2],posc[3]];
                
                V_k1234[c1,spin1,spin2,nb1]=V_k1234[c1,spin1,spin2,nb1]+coe1*coe2*coe3/2;
            end
        end

        
    end
    return V_k1234
end

 
function get_V_k1234_spinful(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234)
    #<k1sigma,k2sigma'|V|k4sigma,k3sigma'>
    Lk=size(k1234,1);
    V_k1234=zeros(Lk,2,2,2)*im;#k1234,sigma,sigma',band
    
    for (spin1, spin2, c1) in Iterators.product(1:2, 1:2, 1:Lk)
        (ind_k1,ind_k2,ind_k3,ind_k4)=k1234[c1,:];
        elem=0;
        for (c2,c3) in Iterators.product(1:2*nK+1, 1:2*nK+1)
            # println(klabel_set[ind_k1,:])
            # println(klabel_set[ind_k4,:])
            dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
            g_label=[c2-1-nK,c3-1-nK];
            q_label=dk_label+g_label;
    
            Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
            # posa=search_label(Fma,q_label,3);
            posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
            if length(posa)>0
                posa=posa[1]
            else
                continue
                #posa=[-1];
            end
            Fmb=Fm_label[ind_k2,ind_k3,:,:,spin2,:];
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

            for nb1 in range(1,2)
                #[posa;posb];
                #posc;
                
                coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb1];
                coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb1,nb1];
                coe3=q_matrix[posc[1],posc[2],posc[3]];
                
                V_k1234[c1,spin1,spin2,nb1]=V_k1234[c1,spin1,spin2,nb1]+coe1*coe2*coe3/2;
            end
        end


    end
    return V_k1234
end


function Hartree_Fock_correction(N,nK,klabel_set,Fm_label,Fm,q_matrix_label,q_matrix)
    
    #define momentum sets for Hatree-fock energies
    k12=Matrix{Int64}(undef,N^2,2)*0;
    step=1;
    for c1=1:N
        for c2=1:N
            k12[step,:]=[c1,c2];
            step=step+1;
        end
    end
    V_Hatree=zeros(N^2,2,2,2,2)*im;#spin1,spin2,band1,band2
    V_Fock=zeros(N^2,2,2,2,2)*im;#spin1,spin2,band1,band2
    #get hatree energies
    for spin1 in range(1,2)
        for spin2 in range(1,2)
            for c1 in range(1,size(k12,1))
                ind_k1=k12[c1,1];
                ind_k2=k12[c1,2];
                ind_k3=k12[c1,2];
                ind_k4=k12[c1,1];
                #k1,k2,k3,k4=k1,k2,k2,k1
                #spin1,spin2,spin3,spin4=spin1,spin2,spin2,spin1
                #nb1,nb2,nb3,nb4=nb1,nb2,nb2,nb1
                
                for c2 in range(1,2*nK+1)
                    for c3 in range(1,2*nK+1)
                        dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
                        g_label=[c2-1-nK,c3-1-nK];
                        q_label=dk_label+g_label;

                        Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
                        #posa=findall(x->x<1e-10,Fma,q_label,3);
                        posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
                        if length(posa)>0
                            posa=posa[1]
                        else
                            continue
                            #posa=[-1];
                        end
                        Fmb=Fm_label[ind_k2,ind_k3,:,:,spin2,:];
                        #posb=findall(x->x<1e-10,Fmb,-q_label,3);
                        posb=findall(x->x<1e-10, abs.(Fmb[:,:,1].+q_label[1])+abs.(Fmb[:,:,2].+q_label[2]))
                        if length(posb)>0
                            posb=posb[1]
                        else
                            continue
                            #posb=[-1];
                        end

                        #posc=findall(x->x<1e-10, q_matrix_label,q_label,4);
                        posc=findall(x->x<1e-10, abs.(q_matrix_label[:,:,:,1].-q_label[1])+abs.(q_matrix_label[:,:,:,2].-q_label[2]))
                        if length(posc)>0
                            posc=posc[1]
                        else
                            continue
                            #posc=[-1];
                        end


                        for nb1=1:2
                            for nb2=1:2        
                                #[posa;posb];
                                #posc;

                                coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb1];
                                coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb2,nb2];
                                coe3=q_matrix[posc[1],posc[2],posc[3]];

                                V_Hatree[c1,spin1,spin2,nb1,nb2]=V_Hatree[c1,spin1,spin2,nb1,nb2]+coe1*coe2*coe3/2;
                            end
                        end

                    end
                end
            end
        end
    end

    #get fock energies
    #with the same spin
    for spin1=1:2
        for c1 in range(1,size(k12,1))
            ind_k1=k12[c1,1];
            ind_k2=k12[c1,2];
            ind_k3=k12[c1,1];
            ind_k4=k12[c1,2];
            #k1,k2,k3,k4=k1,k2,k1,k2
            #spin1,spin2,spin3,spin4=spin1,spin2,spin2,spin1=spin1,spin1,spin1,spin1
            #nb1,nb2,nb3,nb4=nb1,nb2,nb1,nb2

            
            for c2 in range(1,2*nK+1)
                for c3 in range(1,2*nK+1)
                    
                    dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
                    g_label=[c2-1-nK,c3-1-nK];
                    q_label=dk_label+g_label;
                    Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
                    #posa=search_label(Fma,q_label,3);
                    posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
                    if length(posa)>0
                        posa=posa[1]
                    else
                        continue
                        #posa=[-1];
                    end
                    Fmb=Fm_label[ind_k2,ind_k3,:,:,spin1,:];
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

                    
                    for nb1=1:2
                        for nb2=1:2
                            #[posa;posb];
                            #posc;

                            coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb2];
                            coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin1,nb2,nb1];
                            coe3=q_matrix[posc[1],posc[2],posc[3]];

                            V_Fock[c1,spin1,spin1,nb1,nb2]=V_Fock[c1,spin1,spin1,nb1,nb2]+coe1*coe2*coe3/2;
                        end
                    end
                    
                end
            end
        end
    end

    global Hartree_energy_type
    println("Hartree_energy_type= "*Hartree_energy_type);


    #sort Hatree and Fock energies for second band states
    E_Hatree=zeros(N,2)*im;#k,spin
    E_Fock=zeros(N,2)*im;#k,spin
    #V_Hatree(c1,spin1,spin2,nb1=2,nb2=1)
    #V_Fock(c1,spin1,spin1,nb1=2,nb2=1)
    for c1=1:N
        pos1=findall(x->x.==c1,k12[:,1]);
    #     pos2=find(k12(:,2)==c1);
    #     pos=[pos1;pos2];
    #     pos=unique(pos);
        pos=pos1;
        for spin1=1:2
            for cp in range(1,length(pos))
                if Hartree_energy_type=="spinful"
                    E_Hatree[c1,spin1]=E_Hatree[c1,spin1]+V_Hatree[pos[cp],spin1,1,2,1]+V_Hatree[pos[cp],spin1,2,2,1];
                elseif Hartree_energy_type=="spinless"
                    E_Hatree[c1,spin1]=E_Hatree[c1,spin1]+V_Hatree[pos[cp],spin1,spin1,2,1];
                else 
                    error("unknown type");
                end
                #E_Hatree(c1,spin1)=E_Hatree(c1,spin1)+V_Hatree(pos(cp),spin1,spin1,2,1);
                E_Fock[c1,spin1]=E_Fock[c1,spin1]+V_Fock[pos[cp],spin1,spin1,2,1];
            end
        end
        
    end
    E_Hatree=E_Hatree+conj(E_Hatree);
    E_Fock=E_Fock+conj(E_Fock);
    E_Fock=-E_Fock;
    return E_Hatree,E_Fock
end


function Hartree_Fock_correction_model_B(N,nK,klabel_set,Fm_label,Fm,q_matrix_label,q_matrix)
    
    #define momentum sets for Hatree-fock energies
    k12=Matrix{Int64}(undef,N^2,2)*0;
    step=1;
    for c1=1:N
        for c2=1:N
            k12[step,:]=[c1,c2];
            step=step+1;
        end
    end
    V_Hatree=zeros(N^2,2,2,2,2)*im;#spin1,spin2,band1,band2
    V_Fock=zeros(N^2,2,2,2,2)*im;#spin1,spin2,band1,band2
    #get hatree energies
    for spin1 in range(1,2)
        for spin2 in range(1,2)
            for c1 in range(1,size(k12,1))
                ind_k1=k12[c1,1];
                ind_k2=k12[c1,2];
                ind_k3=k12[c1,2];
                ind_k4=k12[c1,1];
                #k1,k2,k3,k4=k1,k2,k2,k1
                #spin1,spin2,spin3,spin4=spin1,spin2,spin2,spin1
                #nb1,nb2,nb3,nb4=nb1,nb2,nb2,nb1
                
                for c2 in range(1,2*nK+1)
                    for c3 in range(1,2*nK+1)
                        dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
                        g_label=[c2-1-nK,c3-1-nK];
                        q_label=dk_label+g_label;

                        Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
                        #posa=findall(x->x<1e-10,Fma,q_label,3);
                        posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
                        if length(posa)>0
                            posa=posa[1]
                        else
                            continue
                            #posa=[-1];
                        end
                        Fmb=Fm_label[ind_k2,ind_k3,:,:,spin2,:];
                        #posb=findall(x->x<1e-10,Fmb,-q_label,3);
                        posb=findall(x->x<1e-10, abs.(Fmb[:,:,1].+q_label[1])+abs.(Fmb[:,:,2].+q_label[2]))
                        if length(posb)>0
                            posb=posb[1]
                        else
                            continue
                            #posb=[-1];
                        end

                        #posc=findall(x->x<1e-10, q_matrix_label,q_label,4);
                        posc=findall(x->x<1e-10, abs.(q_matrix_label[:,:,:,1].-q_label[1])+abs.(q_matrix_label[:,:,:,2].-q_label[2]))
                        if length(posc)>0
                            posc=posc[1]
                        else
                            continue
                            #posc=[-1];
                        end


                        for nb1=1:2
                            for nb2=1:2        
                                #[posa;posb];
                                #posc;

                                coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb1];
                                coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin2,nb2,nb2];
                                coe3=q_matrix[posc[1],posc[2],posc[3]];

                                V_Hatree[c1,spin1,spin2,nb1,nb2]=V_Hatree[c1,spin1,spin2,nb1,nb2]+coe1*coe2*coe3/2;
                            end
                        end

                    end
                end
            end
        end
    end

    #get fock energies
    #with the same spin
    for spin1=1:2
        for c1 in range(1,size(k12,1))
            ind_k1=k12[c1,1];
            ind_k2=k12[c1,2];
            ind_k3=k12[c1,1];
            ind_k4=k12[c1,2];
            #k1,k2,k3,k4=k1,k2,k1,k2
            #spin1,spin2,spin3,spin4=spin1,spin2,spin2,spin1=spin1,spin1,spin1,spin1
            #nb1,nb2,nb3,nb4=nb1,nb2,nb1,nb2

            
            for c2 in range(1,2*nK+1)
                for c3 in range(1,2*nK+1)
                    
                    dk_label=klabel_set[ind_k1,:]-klabel_set[ind_k4,:];
                    g_label=[c2-1-nK,c3-1-nK];
                    q_label=dk_label+g_label;
                    Fma=Fm_label[ind_k1,ind_k4,:,:,spin1,:];
                    #posa=search_label(Fma,q_label,3);
                    posa=findall(x->x<1e-10, abs.(Fma[:,:,1].-q_label[1])+abs.(Fma[:,:,2].-q_label[2]))
                    if length(posa)>0
                        posa=posa[1]
                    else
                        continue
                        #posa=[-1];
                    end
                    Fmb=Fm_label[ind_k2,ind_k3,:,:,spin1,:];
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

                    
                    for nb1=1:2
                        for nb2=1:2
                            #[posa;posb];
                            #posc;

                            coe1=Fm[ind_k1,ind_k4,posa[1],posa[2],spin1,nb1,nb2];
                            coe2=Fm[ind_k2,ind_k3,posb[1],posb[2],spin1,nb2,nb1];
                            coe3=q_matrix[posc[1],posc[2],posc[3]];

                            V_Fock[c1,spin1,spin1,nb1,nb2]=V_Fock[c1,spin1,spin1,nb1,nb2]+coe1*coe2*coe3/2;
                        end
                    end
                    
                end
            end
        end
    end

    global Hartree_energy_type
    println("Hartree_energy_type= "*Hartree_energy_type);


    #sort Hatree and Fock energies for second band states
    E_Hatree=zeros(N,2)*im;#k,spin
    E_Fock=zeros(N,2)*im;#k,spin
    #V_Hatree(c1,spin1,spin2,nb1=2,nb2=1)
    #V_Fock(c1,spin1,spin1,nb1=2,nb2=1)
    for c1=1:N
        pos1=findall(x->x.==c1,k12[:,1]);
    #     pos2=find(k12(:,2)==c1);
    #     pos=[pos1;pos2];
    #     pos=unique(pos);
        pos=pos1;
        for spin1=1:2
            for cp in range(1,length(pos))
                if Hartree_energy_type=="spinful"
                    E_Hatree[c1,spin1]=E_Hatree[c1,spin1]+V_Hatree[pos[cp],spin1,1,2,1]+V_Hatree[pos[cp],spin1,2,2,1];
                elseif Hartree_energy_type=="spinless"
                    E_Hatree[c1,spin1]=E_Hatree[c1,spin1]+V_Hatree[pos[cp],spin1,spin1,2,1];
                else 
                    error("unknown type");
                end

            end
        end
        
    end

    for c1=1:N
        pos1=findall(x->x.==c1,k12[:,2]);
    #     pos2=find(k12(:,2)==c1);
    #     pos=[pos1;pos2];
    #     pos=unique(pos);
        pos=pos1;
        for spin1=1:2
            for cp in range(1,length(pos))
                #E_Hatree(c1,spin1)=E_Hatree(c1,spin1)+V_Hatree(pos(cp),spin1,spin1,2,1);
                E_Fock[c1,spin1]=E_Fock[c1,spin1]+V_Fock[pos[cp],spin1,spin1,2,1];
            end
        end
        
    end
    E_Hatree=E_Hatree+conj(E_Hatree);
    E_Fock=E_Fock+conj(E_Fock);
    E_Fock=-E_Fock;
    return E_Hatree,E_Fock
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