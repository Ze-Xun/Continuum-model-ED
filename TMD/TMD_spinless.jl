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
function main(Nup,Ndn,V1,V2,theta_,flux_1,flux_2,filenm_)

    Ntotal=Nup+Ndn;
    
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
        k2_shift=-flux_2*V1_/N;#check the sign???

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
    println(k_set)
    @assert size(k_set,1)==N;
    klabel_set=k_set;
    # klabel_set=k_set-0.5;



    
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
    V_k1234_general=get_V_k1234_spinup(nK,klabel_set,Fm,Fm_label,q_matrix_label,q_matrix,k1234);
    V_k1234=V_k1234_general[:,:,:,nb];
    #V_k1234=V_k1234_general[:,1,1,nb];
    if simplify_interaction
        println("simplify interaction terms")
        #k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn=reduce_interaction(k1234,V_k1234);
        #k1234_upup,V_k1234_upup=reduce_interaction(k1234,V_k1234);
        k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn=reduce_interaction(k1234,V_k1234);
    end
    ########################################
    



    ##########################################
    basis_up=spinless_fermion_basis_1d(N,Nup);
    ########################################



    n_E=10;
    occu_set_up,particle_set_up=get_spinless_occu(N,Nup,basis_up,length(basis_up));

    k_total_set=klabel_set;
    # k_total_set[:,1]=k_total_set[:,1].-k1_shift*(Nup+Ndn-1);
    # k_total_set[:,2]=k_total_set[:,2].-k2_shift*(Nup+Ndn-1);
    ksector_size=zeros(Int64,size(k_total_set,1));

    @tensor Kset_up[:]:=occu_set_up[-1,1]*klabel_set[1,-2];


    ########################
    println("start ED:");flush(stdout);

    #Qn_set=list();
    #Qn_set.append(list([Nup,Ndn]));
    En_set=Vector{Vector}(undef,N);


    for nsec in range(1,N)
        starting_time=now()
        ksector_size,ksector_ind=get_sector_dim_spinup(N,basis_up,Kset_up,k_total_set,nsec);
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
        spin1=1;
        global break_Hamiltonian;
        Hk=kinetic_energy_1band_spinup(N,ksector_ind,E_band[:,spin1,nb]);
        global use_kinetic_energy;
        println("use_kinetic_energy="*string(use_kinetic_energy));
        if use_kinetic_energy==false
            Hk=Hk*0;
        end

        
        if ~simplify_interaction
            # V_upup=interaction_1band_spinup(N,ksector_ind, k1234,V_k1234);
        else
            if ~break_Hamiltonian
                #k1234_upup,V_k1234_upup, k1234_dndn,V_k1234_dndn,  k1234_updn,V_k1234_updn
                # V_upup=interaction_1band_spinup(N,ksector_ind, k1234,V_k1234);
            else
                H_set=Vector{SparseMatrixCSC}(undef,0);
                #push!(H_set,Hk);
                push!(H_set,ComplexF32.(Hk));
            
                H_set=interaction_1band_spinup_break_H(H_set,N,ksector_ind, k1234_upup,V_k1234_upup);
            end
        end
 
        if Dim>2000000
            GC.gc(true);
        end
        if ~break_Hamiltonian
            H=V_upup+Hk;
        else
            #H=Hk+V_upup;
            function apply_H(x::Vector,H_set::Vector{SparseMatrixCSC})
                x1=H_set[1]*x;
                for cc=2:length(H_set)
                    x1=x1+H_set[cc]*x;
                end
                return x1
            end
            F_H(x)=apply_H(x,H_set);
            
        end

        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for construct Hamiltonian: "*string(Time));flush(stdout);

        starting_time=now()
        #eu,ev=H.eigh()
        if ~break_Hamiltonian
            if Dim<10
                eu,ev=eigen(Matrix(H));
            else
                #println("compare:")
                #eu,ev=eigs(H,nev=min(n_E,Dim-2),which=:SR)
                eu,ev=eigsolve(H,min(n_E,Dim-2),:SR, ishermitian=true)
            end
        else
            GC.gc(true);
            x_init=ComplexF32.(randn(Dim)*im);x_init=x_init/norm(x_init);
            #x_init=randn(Dim)*im;x_init=x_init/norm(x_init);
            println("test speed of Mv");flush(stdout);
            @time F_H(x_init);
            @time F_H(x_init);
            eu,ev=eigsolve(F_H,x_init,min(n_E,Dim-2),:SR, ishermitian=true)
        end
        @assert norm(imag.(eu))/norm(eu)<1e-12
        eu=real.(eu);
        eu=sort(eu);
        println(eu)
        
        Now=now();
        Time=Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(Now) - Dates.DateTime(starting_time)));
        println("time for diagonalizing: "*string(Time));flush(stdout);
        
        En_set[nsec]=eu;
        GC.gc(true);
    end


    matwrite(filenm_, Dict(
        "k_total_set" => k_total_set,
        "En_set" => En_set
    ); compress = false)

end
# ##########################
band_model=band_single_k;#band_single_k
global use_kinetic_energy,break_Hamiltonian,band_model,simplify_interaction
use_kinetic_energy=true;
break_Hamiltonian=true;#break big Hamiltonian into sections so that construction sparse Hamiltonian is possible
simplify_interaction=true;
Nup=2;
Ndn=0;
theta_=2;
flux_1=0;
flux_2=0;
filenm_ = "FQHE_"*"theta"*string(theta_)*"_Nup"*string(Nup)*"Ndn"*string(Ndn)*"_flx_"*string(flux_1)*"_"*string(flux_2)*".mat";

#cluster vector
# V1=[6,-3];#integer vector
# V2=[3,-6];#integer vector
V1=[4,0];#integer vector
V2=[0,4];#integer vector

main(Nup,Ndn,V1,V2,theta_, flux_1,flux_2, filenm_)

# flux_length=10;
# for cc = 1:flux_length
#     simplify_interaction=true;
#     Nup=5;
#     Ndn=0;
#     theta_=2;
#     flux_1=cc/flux_length;
#     flux_2=0;
#     filenm_ = "FQHE_"*"theta"*string(theta_)*"_Nup"*string(Nup)*"Ndn"*string(Ndn)*"_flx_"*string(flux_1)*"_"*string(flux_2)*".mat";
#     main(Nup,Ndn,theta_, flux_1,flux_2, simplify_interaction, filenm_)
# end


