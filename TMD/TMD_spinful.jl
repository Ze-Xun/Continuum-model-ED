using Revise, TensorKit
using JLD2,MAT
using KrylovKit
using JSON
using Random
using Dates
using SparseArrays,Combinatorics
#using Arpack
using KrylovKit

cd(@__DIR__)
include("ED_lib.jl")
include("model_lib.jl")

####################
import LinearAlgebra.BLAS as BLAS
n_cpu=4;
BLAS.set_num_threads(n_cpu);
println("number of cpus: "*string(BLAS.get_num_threads()))
Base.Sys.set_process_title("C"*string(n_cpu)*"_ED")
pid=getpid();
println("pid="*string(pid));
####################

function main(band_model,Nup,Ndn,theta_,flux_type,flux_1,flux_2,simplify_interaction,filenm_)




    Ntotal=Nup+Ndn;
    #@assert (Nup>0)&(Ndn>0);
    #cluster vector
    # V1=[6,-3];#integer vector
    # V2=[3,-6];#integer vector
    V1=[4,0];#integer vector
    V2=[0,4];#integer vector

    @assert 0<=flux_1<=1
    @assert 0<=flux_2<=1
    if (abs(flux_1)>0)|(abs(flux_2)>0)
        @assert (V1[1]==0)|(V1[2]==0);
        @assert (V2[1]==0)|(V2[2]==0);
    end


    N_bands=2;
    spin=1;#1 or 2
    nK=6;#number of g vectors in each side


    H,eu,ev,g_set,kapap,kapam,basis,H_x,H_y,g_lattice=band_single_k(theta_,([0,0],[0,0]),nK);

    #####################################

    g1=g_set[1,:];
    g3=g_set[3,:];
    #g_lattice: multiples of g vector

    g1_3d=[g1[1],g1[2], 0];
    g3_3d=[g3[1],g3[2], 0];
    z_3d=[0,0,1];

    a1=2*pi*cross(z_3d,g3_3d)/(dot(g1_3d, cross(z_3d,g3_3d)));
    a2=2*pi*cross(z_3d,g1_3d)/(dot(g3_3d, cross(z_3d,g1_3d)));


    V1_3d=V1[1]*a1+V1[2]*a2;
    V2_3d=V2[1]*a1+V2[2]*a2;

    N_=abs(cross(V1_3d,V2_3d)[3]/cross(a1,a2)[3]);
    @assert (Int(abs(round(N_)))-abs(N_))/abs(round(N_))<1e-14;
    N=Int(abs(round(N_)));
    println("cluster size N="*string(N))
    G1_3d=2*pi*cross(z_3d,V2_3d)/(dot(V1_3d, cross(z_3d,V2_3d)));
    G2_3d=2*pi*cross(z_3d,V1_3d)/(dot(V2_3d, cross(z_3d,V1_3d)));

    G1=G1_3d[1:2];
    G2=G2_3d[1:2];


    #################

    if (V1[2]==0)&&(V2[1]==0)
        k_set=zeros(N,2);
        Step=0;
        for c2=1:V2[2]
            for c1=1:V1[1]
                label=[(c1-1)/V1[1],(c2-1)/V2[2]];
                Step=Step+1;
                k_set[Step,:]=label;
            end
        end

        #determine shift of k points
        #old method
        # k1_shift=flux_1*[1/abs(V1[1]),0];
        # k2_shift=flux_2*[0,1/abs(V2[2])];

        #new method
        V1_=[V1[2],-V1[1]];
        V2_=[V2[2],-V2[1]];
        k1_shift=flux_1*V2_/N;
        k2_shift=-flux_2*V1_/N; #check the sign???
    else

        LL=40;
        k_set=zeros((LL*LL,2));
        Step=0;
        for c2 in range(1,LL)
            for c1 in range(1,LL)
                p1=Int(c1-round(LL/2));
                p2=Int(c2-round(LL/2));

                V1_=[V1[2],-V1[1]];
                V2_=[V2[2],-V2[1]];
                label=[(p1*V2_[1]+p2*V1_[1])//N, (p1*V2_[2]+p2*V1_[2])//N];
                
                if (0<=label[1])&(label[1]<1)&(0<=label[2])&(label[2]<1)
                    Step=Step+1;
                    k_set[Step,:]=label;
                end
            end
        end
        k_set=k_set[1:Step,:];
        #determine shift of k points
        k1_shift=flux_1*V2_/N;
        k2_shift=-flux_2*V1_/N;#check the sign???
    end
    # println(k_set)
    @assert size(k_set,1)==N;
    klabel_set=k_set;
    # klabel_set=k_set-0.5;




    # klabel_set[:,1]=klabel_set[:,1].-k1_shift;
    # klabel_set[:,2]=klabel_set[:,2].-k2_shift;


    ########################
    momentum_set,band_state,E_band=get_bandstate(band_model,N,nK,N_bands,g1,g3,theta_,klabel_set,k1_shift,k2_shift,"spin_independent")
    println("k points:")
    println(momentum_set)
    println("E_band:")
    println(E_band[:,:,1]);flush(stdout);
    Fm, Fm_label=get_Fm(N,nK,N_bands,band_state,momentum_set,klabel_set);



    #coulomb interaction
    q_matrix=zeros((N,2*nK+1,2*nK+1));
    q_matrix_label=zeros((N,2*nK+1,2*nK+1,2));
    a0=3.52e-10;
    theta=theta_/180*pi;
    aM=a0/(2*sin(theta/2));
    k0 = 8.98755e9;
    J_to_meV = 6.24150636e21; #meV
    e_charge = 1.60217663e-19; #coulomb
    epsilon_r = 10.0;
    Area = sqrt(3)/2*N*aM^2;
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
                    q_matrix[c1,c2,c3] = 2*pi*e_charge^2/(epsilon_r*sqrt(qQ[1]^2+qQ[2]^2))*J_to_meV/Area*k0;
                end
            end
        end
    end

        
    #################################



    #cdag_{k1,sigma} cdag_{k2,sigma'} c_{k3,sigma'} c_{k4,sigma}
    Step=0;
    k1234=zeros(Int,(N^4,4));
    for c1 in range(1,N)
        for c2 in range(1,N)
            for c3 in range(1,N)
                for c4 in range(1,N)
                    k_label_sum=sum(klabel_set[[c1,c2],:]-klabel_set[[c3,c4],:],dims=1);
                    k_label_sum_int=round.(k_label_sum);
                    if (abs(k_label_sum[1]-k_label_sum_int[1])<1e-10)&(abs(k_label_sum[2]-k_label_sum_int[2])<1e-10)
                        Step=Step+1;
                        k1234[Step,:]=[c1,c2,c3,c4];
                    end
                end
            end
        end
    end
                        
                        
    k1234=k1234[1:Step,:];

    ########################################
    nb=1;#index of band
    V_k1234_general=get_V_k1234_spinful(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234);
    V_k1234=V_k1234_general[:,:,:,nb];
    if simplify_interaction
        println("simplify interaction terms")
        k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn=reduce_interaction(k1234,V_k1234);
    end
    ########################################


    ##########################################
    basis_up=spinless_fermion_basis_1d(N,Nup);
    basis_dn=spinless_fermion_basis_1d(N,Ndn);
    ##########################################



    n_E=10;
    occu_set_up,particle_set_up=get_spinless_occu(N,Nup,basis_up,length(basis_up));
    occu_set_dn,particle_set_dn=get_spinless_occu(N,Ndn,basis_dn,length(basis_dn));


    k_total_set=klabel_set;
    # k_total_set[:,1]=k_total_set[:,1].-k1_shift*(Nup+Ndn-1);
    # k_total_set[:,2]=k_total_set[:,2].-k2_shift*(Nup+Ndn-1);
    ksector_size=zeros(Int64,size(k_total_set,1));



    @tensor Kset_up[:]:= occu_set_up[-1,1]*klabel_set[1,-2];
    @tensor Kset_dn[:]:= occu_set_dn[-1,1]*klabel_set[1,-2];


    ########################
    println("start ED:");flush(stdout);

    #Qn_set=list();
    #Qn_set.append(list([Nup,Ndn]));
    En_set=Vector{Vector}(undef,N);


    #############################################





    for nsec in range(1,N)

        starting_time=now()
        ksector_size, ksector_total_ind, ksector_ind=get_sector_dim_spinful(N,basis_up,basis_dn,Kset_up,Kset_dn,k_total_set,nsec);

        Dim=ksector_size;
        println("dim="*string(Dim));flush(stdout);
        # ksector_ind=ksector_ind[end:-1:1]
        # println(ksector_size)
        # println(ksector_total_ind)
        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for construct basis: "*string(Time));flush(stdout);
        starting_time=now()
        #basis,H=ED_spinup(N,Nup,external_basis,H_terms)
        
        Hk=kinetic_energy_1band_spinful(N,ksector_total_ind,E_band[:,:,nb]);

        nthreads = Threads.nthreads()
        println("number of threads: "*string(nthreads));flush(stdout);
        
        if ~simplify_interaction
            spin1=1;
            spin2=1;
            V_upup=interaction_1band_spinful(N,ksector_total_ind, k1234,V_k1234,spin1,spin2);
            spin1=1;
            spin2=2;
            V_updn=interaction_1band_spinful(N,ksector_total_ind, k1234,V_k1234,spin1,spin2);
            spin1=2;
            spin2=1;
            V_dnup=interaction_1band_spinful(N,ksector_total_ind, k1234,V_k1234,spin1,spin2);
            spin1=2;
            spin2=2;
            V_dndn=interaction_1band_spinful(N,ksector_total_ind, k1234,V_k1234,spin1,spin2);

            if Nup==0
                H=Hk+V_dndn;
            elseif Ndn==0
                H=Hk+V_upup;
            else
                H=Hk+V_upup+V_updn+V_dnup+V_dndn;
            end
        else
            #k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn
            
            spin1=1;
            spin2=1;
            V_upup=interaction_1band_spinful(N,ksector_total_ind, k1234_upup,V_k1234_upup,spin1,spin2);

            spin1=1;
            spin2=2;
            V_updn=interaction_1band_spinful(N,ksector_total_ind, k1234_updn,V_k1234_updn,spin1,spin2);

            spin1=2;
            spin2=2;
            V_dndn=interaction_1band_spinful(N,ksector_total_ind, k1234_dndn,V_k1234_dndn,spin1,spin2);

            if Nup==0
                H=Hk+V_dndn;
            elseif Ndn==0
                H=Hk+V_upup;
            else
                H=Hk;
                if maximum(size(V_upup))>0
                    H=H+V_upup;
                end
                if maximum(size(V_dndn))>0
                    H=H+V_dndn;
                end
                if maximum(size(V_updn))>0
                    H=H+V_updn;
                end
            end
        end



        

        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for construct Hamiltonian: "*string(Time));flush(stdout);

        starting_time=now()
        #eu,ev=H.eigh()
        
        if Dim<10
            eu,ev=eigen(Matrix(H));
        else
            #println("compare:")
            #eu,ev=eigs(H,nev=min(n_E,Dim-2),which=:SR)
            eu,ev=eigsolve(H,min(n_E,Dim-2),:SR, ishermitian=true)
        end
        @assert norm(imag.(eu))/norm(eu)<1e-12
        eu=real.(eu);
        eu=sort(eu);
        println(eu)
        
        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for diagonalizing: "*string(Time));flush(stdout);
        
        En_set[nsec]=eu;



    end

    matwrite(filenm_, Dict(
        "k_total_set" => k_total_set,
        "En_set" => En_set
    ); compress = false)
end
##########################



# ##########################
band_model=band_single_k;#band_single_k
simplify_interaction=true;
Nup=2;
Ndn=2;
theta_=2;
flux_1=0;
flux_2=0;
flux_type="spin_dependent";#"spin_dependent" or "spin_independent" 
filenm_ = "FQHE_"*"theta"*string(theta_)*"_Nup"*string(Nup)*"Ndn"*string(Ndn)*"_flx_"*string(flux_1)*"_"*string(flux_2)*".mat";
main(band_model,Nup,Ndn,theta_, flux_type, flux_1,flux_2,  simplify_interaction, filenm_)

# flux_length=10;
# for cc in range(1,flux_length)
#     simplify_interaction=true;
#     Nup=5;
#     Ndn=0;
#     theta_=2;
#     flux_1=cc/flux_length;
#     flux_2=0;
#     filenm_ = "FQHE_"*"theta"*string(theta_)*"_Nup"*string(Nup)*"Ndn"*string(Ndn)*"_flx_"*string(flux_1)*"_"*string(flux_2)*".mat";
#     main(Nup,Ndn,theta_, flux_1,flux_2, flux_type, simplify_interaction, filenm_)
# end
