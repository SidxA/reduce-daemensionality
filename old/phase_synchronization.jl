#==================================================================#
#               PLAN: THEORY
#==================================================================#
#data: big vector                                   needs specified matrix          T x N
#   centralize in place
#   PCA by svd
#eigenvectors: choose a number                      needs specified matrix          T x K
#   selection
#   again(?) centralize
#   SSA(global) select W                            needs specified lag-matrix      K x W x P(=T-W+1)
#modes : here we only want the season               needs specified matrix          K x T x 2
#   deduce the season in place
#   SSA(global) select L, on the rest               needs specified lag-matrix      K x L x Q(=T-L+1)
#modes                                              needs specified matrix          K x T x M
#   assign phase variable to each                   needs specified matrix          K x T x M

# variables
# T series length
# N number of original variables
# K chosen number of kept core variables
# W chosen window length of seasonalizing SSA
# P number of seasonalizing SSA windows
# L chosen window length of descaling SSA
# Q number of descaling SSA windows
# M chosen number of different time scale modes

#==================================================================#
#               PLAN: EXECUTION
#==================================================================#
# outer class
#   initialize with the numbers T,N,K,W,P,L,Q,M
#
#   constructor method to create the specific matrices
#
#   outer method that calls the individual methods
#       centralize data individually
#       pcaand take first couple of pc for reconstruction
#       iterate over the selected, centralized eigenvectors
#           perform SSA for season
#           take the season from the eigenvector
#           perform SSA for scale
#           iterate over the chosen number of scale bands
#               estimate the phase maybe use DMD            !!! hilbetrt trafo empiric wavemode decomposition
#
#   individual: PCA             take data and pca matrix and write to the pca matrix
#   individual: global SSA      take one row of pca matrix and delay-matrix and write to delay-matrix
#   individual: centralize      take any row and write to this row
#   individual: phase-assign

#        REQUIRE THE FIRST LAST YEAR PROLONGATION AND LATER CUTOFF
#        THINK OF WHERE TO PUT IN THE CENTRALIZE

#==================================================================#
#               class test
#==================================================================#
#requirements
using Statistics
using Setfield # maybe change to immutable later
using MultivariateStats
using Random
using CSV
using DataFrames


mutable struct algo
    #------------------------------------------------------------------
    # variables
    #------------------------------------------------------------------
    T::Int# T series length
    N::Int# N number of original variables
    N_sep::Int # separation mark : atmosphere up to this point
    K::Int# K chosen number of kept core variables from BOTH series
    W::Int# W chosen window length of seasonalizing SSA
    P::Int# P number of seasonalizing SSA windows
    L::Int# L chosen window length of descaling SSA
    Q::Int# Q number of descaling SSA windows
    M::Int# M chosen number of different time scale modes
    data            #T x N
    eigenvectors    #T x K
    eigvec_preseas #T x K
    delay_trj_seas  #K x W x P
    modes_seas      #K x T x 2
    delay_trj_scal  #K x L x Q
    modes_scal      #K x T x M
    phases          #K x T x M

    pca        #placeholder
    seas       #placeholder W x P
    ph1        #placeholder
    ph2        #placeholder

    #------------------------------------------------------------------
    # initialize by prior knowledge of the size of the data
    #------------------------------------------------------------------
    function algo(T,N,N_sep)

        K = 6
        W = Int(floor(T/3))
        P = T - W + 1
        L = Int(floor(T/3))
        Q = T - L +1
        M = 6
        data            = Array{Float64}(undef,T,N)
        eigenvectors    = Array{Float64}(undef,T,K)
        eigvec_preseas = Array{Float64}(undef,T,K)
        delay_trj_seas  = Array{Float64}(undef,K,W,P)
        modes_seas      = Array{Float64}(undef,K,T,2)
        delay_trj_scal  = Array{Float64}(undef,K,L,Q)
        modes_scal      = Array{Float64}(undef,K,T,M)
        phases          = Array{Float64}(undef,K,T,M)
        seas            = Array{Float64}(undef,W,P)

        return new(T,N,N_sep,K,W,P,L,Q,M,
                    data,eigenvectors,eigvec_preseas,delay_trj_seas,modes_seas,
                    delay_trj_scal,modes_scal,phases,PCA,seas,:None,:None)
    end

end



    #------------------------------------------------------------------
    #helpers
    #------------------------------------------------------------------

            #------------------------------------------------------------------
            #centralize
            #------------------------------------------------------------------
            function centralize(data)
                m = mean(data)
                cv = std(data)
                return  (data.-m)./cv
            end
            #------------------------------------------------------------------
            #load the data                  DUMMY THIS NEEDS TO KNOW THE FILTERED RESULT PRIOR IN THE INIT STEP
            #------------------------------------------------------------------
            function load_the_data(class,name,sep)

                df = DataFrame(CSV.File(name,delim=",",header=2))
                #remove timestamp
                df = df[:,2:end]
                #remove constant rows
                il=[]
                for i in 1:ncol(df)
                    if df[2,i]==df[3,i]
                        il = push!(il,false)
                    elseif any(i -> abs(i) > 500,df[:,i])
                        il = push!(il,false)
                    else
                        il = push!(il,true)
                    end
                end

                df = Array(df)
                il = Array(Bool.(il))
                df = df[:,il]

                # NEED TO ADD FIRST AND LAST YEAR
                
                setfield!(class,:T,size(df)[1])
                setfield!(class,:N,size(df)[2])
                setfield!(class,:N_sep,sep)
                setfield!(class,:data,df)

                return  nothing
            end
            #------------------------------------------------------------------
            #eigenvectors by PCA                UGLY NUMBER HACK FOR ATMO/BIO; NEEDS BETTER DATA STRUCTURE AS NOT TO TRANSPOSE ALL THE TIME #savepoint might need to be inside
            #------------------------------------------------------------------
            function eigvec_by_PCA(class)
                setfield!(class,:pca,fit(PCA,centralize(class.data[:,1:class.N_sep]'),method=:cov,maxoutdim=Int(class.K/2)))
                class.ph1 = MultivariateStats.transform(class.pca,class.data[:,1:class.N_sep]')
                setfield!(class,:pca,fit(PCA,centralize(class.data[:,class.N_sep+1:end]'),method=:cov,maxoutdim=Int(class.K/2)))
                class.ph2 = MultivariateStats.transform(class.pca,class.data[:,class.N_sep+1:end]')
                setfield!(class,:eigenvectors,vcat(class.ph1,class.ph2)')

                return nothing
            end

            #------------------------------------------------------------------
            #de_seasonalizing by global SSA     UGLY NUMBER HACK FOR THE COLUMN OF THE EIGENVECTORS MATRIX
            #------------------------------------------------------------------
            function de_seasonalizing_global_SSA(class)
                setfield!(class,:pca,:None)
                setfield!(class,:eigvec_preseas,class.eigenvectors)
                for i in 1:class.K
                    setfield!(class,:seas,embed_lag(class,i,class.W))
                    setfield!(class,:pca,fit(PCA,class.seas,method=:cov,maxoutdim=2))

                    #class.modes_seas[i,:,:] = transform_global_SSA(class,i,class.W)
                    class.modes_seas[i,:,:] = transform_global_SSA(i,class.W,class.pca,class.eigenvectors,class.T,class.P,class.seas)

                    class.eigenvectors[:,i] = centralize(class.eigenvectors[:,i]) .- centralize(sum(class.modes_seas[i,:,:],dims=2))

                end

                #@set class.delay_trj_seas[number,:,:] = embed_lag(class,number,class.W)
                #@set class.pca = fit(PCA,class.delay_trj_seas[number,:,:],method=:cov,maxoutdim=2)

                #@set class.modes_seas[number,:,:] = transform_global_SSA(class,number,class.W)
                #@set class.eigenvectors[:,number] = class.eigenvectors[:,number] .- sum(class.modes_seas,dims=3)

                return  nothing
            end

            #------------------------------------------------------------------
            #de_scaling by global SSA
            #------------------------------------------------------------------
            function de_scaling_global_SSA(class)
                setfield!(class,:pca,:None)
                for i in 1:class.K
                    setfield!(class,:seas,embed_lag(class,i,class.L))
                    setfield!(class,:pca,fit(PCA,class.seas,method=:cov,maxoutdim=class.M))

                    #class.modes_seas[i,:,:] = transform_global_SSA(class,i,class.W)
                    class.modes_scal[i,:,:] = transform_global_SSA(i,class.L,class.pca,class.eigenvectors,class.T,class.Q,class.seas)

                end

                #@set class.delay_trj_seas[number,:,:] = embed_lag(class,number,class.W)
                #@set class.pca = fit(PCA,class.delay_trj_seas[number,:,:],method=:cov,maxoutdim=2)

                #@set class.modes_seas[number,:,:] = transform_global_SSA(class,number,class.W)
                #@set class.eigenvectors[:,number] = class.eigenvectors[:,number] .- sum(class.modes_seas,dims=3)

                return  nothing
            end

            #------------------------------------------------------------------
            #embedding for global SSA            UGLY BECAUSE IT STILL NEED MANUAL CHOICE BETWEEN THE WINDOW LENGHTES   THIS IS SOME DynAMICAL CENTRALIZING
            #------------------------------------------------------------------
            function embed_lag(class,number,W)
                Y = []
                for i=1:class.T-W+1
                    Y = append!(Y,centralize(class.eigenvectors[i:i+W-1,number]))
                end
                return reshape(float.(Y),W,class.P) # K W P 
            end

            #------------------------------------------------------------------
            #reconstruct components from global SSA   UGLY BECAUSE IT STILL NEED MANUAL CHOICE BETWEEN THE WINDOW LENGHTES   
            #------------------------------------------------------------------
            function transform_global_SSA(number,W,pca,eigenvectors,T,P,data)
                println(number)
                # for both seasons!
                
                for season = 1:2
                    EOF = vec(projection(pca)[:,season])#we have the EOF of length W
                    PC = vec(MultivariateStats.transform(pca,data)[season,:])#we have the PC of length P=N-W+1
                    #we have the lagged trajectory of size (W,P)
                
                    N = T
                    P = P
                
                    RC = []
                    RC_l = []
                    M,L,U = 0.0, 0.0, 0.0
                
                    #first part
                    for t=1:W-1
                        M = 1/t 
                        L = 1
                        U = t
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                
                    #second part
                    for t=W:P
                        M = 1/W
                        L = 1
                        U = W
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                
                    #third part
                    for t=P+1:N
                
                        M = 1/(N-t+1)
                        L = t-N+W
                        U = W
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                
                    #@set class.modes_seas[number,:,season] = RC
                    RC_l = append!(RC_l,RC)
                end
                return reshape(RC_l,T,2) # T,2
            end

            #------------------------------------------------------------------
            #reconstruct components for rescaling for global SSA   UGLY BECAUSE IT STILL NEED MANUAL CHOICE BETWEEN THE WINDOW LENGHTES   


            # UUUUGLY
            
            #------------------------------------------------------------------
            function transform_global_SSA_for_scale(number,L,pca,eigenvectors,T,Q,data)
                println(number)
                # for both seasons!
                RC_l = []
                    EOF = vec(projection(pca)[:,mode])#we have the EOF of length L
                    PC = vec(MultivariateStats.transform(pca,data)[season,:])#we have the PC of length P=N-W+1
                    #we have the lagged trajectory of size (W,P)
                
                    N = T
                    P = P
                
                    RC = []
                    RC_l = []

                    M,L,U = 0.0, 0.0, 0.0
                
                    #first part
                    for t=1:W-1
                        M = 1/t 
                        L = 1
                        U = t
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                
                    #second part
                    for t=W:P
                        M = 1/W
                        L = 1
                        U = W
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                
                    #third part
                    for t=P+1:N
                
                        M = 1/(N-t+1)
                        L = t-N+W
                        U = W
                        inner_sum=[]
                        for j in L:U 
                            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
                        end
                        RC = push!(RC,M*sum(inner_sum))
                    end
                # took away an "end" here that might have not been quite well placed
            
                    #@set class.modes_seas[number,:,season] = RC
                    RC_l = append!(RC_l,RC)
                
                return reshape(RC_l,T,2) # T,2
            end

#==================================================================#
# function to extract the variable
#==================================================================#

# ISOMAP approach 1: fixed number/radius
#   connectivity graph
#   shortest interpoint distances
#   minimize deviation from this to the SVD of theoriginal metric space
#   of the matrix of inner products tau of each

#crude: fit ISOMAP from ManifoldLearning.jl, parameter k


# ISOMAP approach 2: precompute for different thresholds, then choose this
# (online) PCA aka Guido approach

#==================================================================#
# function to decompose into the modes
#==================================================================#

# ssa approach 1: averaging
#   the tough K=N-W+1 all slices approach and a window length of about half

# implement a crude delay window function

function delay_trajectory(X,W)
    N = length(X)
    Y = []
    for i=1:N-W+1
        Y = append!(Y,X[i:i+W-1])
    end
    return reshape(float.(Y),W,N-W+1)
end

# perform the SSA globally on this window
# best make a class global_ssa

function global_SSA(X)
    #this needs to give a matrix of the eigenvectors EOF and a vector of their variance
    #this needs to have a specified size
    # best would be to insert this in a function that assigns a space and then fills it with these numbers
    # possibly use the advanced pca for this from guido
    pca = fit(PCA,X,method=:auto,maxoutdim=50)
    return principalvars(pca),transform(pca,X)
end

function reconstruct(T,W,EOF,PC)

    #we have the EOF of length W
    #we have the PC of length P=N-W+1
    #we have the lagged trajectory of size (W,P)

    N = T
    P = N-W+1

    RC = []
    M,L,U = 0.0, 0.0, 0.0

    #first part
    for t=1:W-1
        M = 1/t 
        L = 1
        U = t
        inner_sum=[]
        for j in L:U 
            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
        end
        RC = push!(RC,M*sum(inner_sum))
    end

    #second part
    for t=W:P
        M = 1/W
        L = 1
        U = W
        inner_sum=[]
        for j in L:U 
            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
        end
        RC = push!(RC,M*sum(inner_sum))
    end

    #third part
    for t=P+1:N

        M = 1/(N-t+1)
        L = t-N+W
        U = W
        inner_sum=[]
        for j in L:U 
            inner_sum = append!(inner_sum, PC[t-j+1]*EOF[j])
        end
        RC = push!(RC,M*sum(inner_sum))
    end

    return RC
end

# ssa approach 2: localized repetitions
#   repeat first and last year
#   moving window of about a year
#   ssa in each window seperately by repeating 5 times
#   identify the annual–seasonal SSA component, ignore the rest
#   reconstruction by weighted sum: each point is covered by 2 windows and they are weighted by the sin/bell
#   merge to the complete time series for annual-seasonal
#   deseasonalize the series using these modes
#   3 years embedding window for low frequency modes with standard SSA
#   cut at 4.5 month, rest is high frequency

#==================================================================#
# function to sort the modes, separate the timescales
#==================================================================#



#==================================================================#
# function to test for phase synchrony
#==================================================================#

#   assign phase variables by local differentiation or by fitted sine




#time stamp for saving analysis results
using Dates


function init_logging()
    weekno = week(unix2datetime(time()))
    datestring = string("KW_",lpad(weekno,2,"0"),"/")
    workdir = "/net/home/lschulz/logs/"
    dir = workdir*datestring
    if isdir(dir)==false
        mkdir(dir)
    end
    return dir
end

dir = init_logging()

#dummy training run

a = algo(7304,92,20)
using Plots
load_the_data(a,"testdata/FLX_US-PFa_FLUXNET2015_FULLSET_DD_1995-2014_1-4.csv",20)
#savefig(heatmap(Matrix(a.data)),dir*"data")
#eigvec_by_PCA(a)
#savefig(plot(a.eigenvectors),dir*"eigenvectors")
#for i in 1:a.K
#    a.eigenvectors[:,i] = centralize(a.eigenvectors[:,i])
#end
#de_seasonalizing_global_SSA(a)
#savefig(plot(a.eigenvectors),dir*"eigenvectors_centralized")
#savefig(plot(a.eigenvectors),dir*"eigenvectors_de_seas")
#savefig(plot([a.eigenvectors[:,1] a.eigvec_preseas[:,1]]),dir*"comparison")
#de_scaling_global_SSA(a)

dir2 = "/net/home/lschulz/logs/KW_18/eemdrun/"
mkdir(dir2)
#mkdir(dir2)
x = centralize(a.data[:,20])
for n = Int.(2:2:8), Sifting_number = Int.(2:2:8), n_IMF = 10, S = Int.(2:2:8)
    title_x = "n "*string(n)*" sift "*string(Sifting_number)*" n_IMF "*string(n_IMF)*" S "*string(S)
    x_emd = eemd(x,EEMDSetting(n,Sifting_number,S,n_IMF))
    savefig(heatmap([x x_emd],title=title_x),dir2*"heat_"*title_x)
end