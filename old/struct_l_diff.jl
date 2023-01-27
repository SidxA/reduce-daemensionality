mutable struct local_diffusion
    #------------------------------------------------------------------
    # variables
    #------------------------------------------------------------------
    T::Int# T series length

    rep_N::Int # number of repetitions
    loc_rate::Float64 # localization rate
    loc_win_len::Int # localization length

    loc_len::Int #len with repetitions
    loc_N::Int #number of localization windows

    W::Int# W chosen window length
    P::Int# P number of windows

    k::Int# k chosen modes
    a::Float64
    e::Int
    t::Int
    #------------------------------------------------------------------
    # placeholders
    #------------------------------------------------------------------
    data            #T
    localized_data  #loc_N x loc_len

    kept_RC         # loc_N x loc_len x k

    #places inside the localize loop where we only keep the RC at the end
    EOF     # W x k
    PC      # k x P
    loc_emb         # P x W
    #------------------------------------------------------------------
    # initialize by putting in the localized data and parameters, maybe just redo with this?
    #------------------------------------------------------------------
    #data, loc_win_len, loc_rate, rep_N, W, k,a,e,t
    function local_diffusion(
        signal::Vector,
        loc_win_len,
        loc_rate,
        rep_N,
        W,
        k,
        a,
        e,
        t
        )

        @assert isodd(rep_N) "replication not odd"

        #embedding part

        data = isodd(length(signal)) ? signal[1:end-1] : signal
        T = length(signal)

        loc_len = loc_win_len * rep_N                               #len with repetitions
        P = loc_len - W + 1                                         # P number of windows

        @assert ((loc_len+1)/2>W) "choose smaller W"

        loc_N = floor(Int,((T+1)/(loc_win_len * loc_rate)))
        localized_data = zeros(loc_N,loc_len)

            #println("replicate loc_emb of signal_l $T")
            #println("windows $loc_size, cover_rate $loc_rate, N_repl $N_repl, throwing away $(Int(T-N_locs*loc_size*loc_rate))")
            ind = 1
            for t in 1:loc_N
                ind_l = ind + loc_win_len - 1 #position_end is position + loc_win_len -1
                if ind_l < T #if position_end in the series
                    localized_data[t,:] = replicate(signal[ind:ind_l],rep_N) #put the window and put in the replications
                    ind += Int(loc_win_len * loc_rate) #position shifts by loc_win_size / loc_rate
                else
                    break
                end
            end


        #everything else
        kept_RC         = Array{Float64}(undef,loc_N,loc_len,k)     # loc_N x loc_len x k
        loc_emb         = Array{Float64}(undef,W,P)                 # W x P
        EOF             = Array{Float64}(undef,W,k)                 # W x k
        PC              = Array{Float64}(undef,k,P)                 # k x P

        return new(
            T,# T series length
            rep_N, # number of repetitions
            loc_rate, # localization rate
            loc_win_len, # localization length
            loc_len, #len with repetitions
            loc_N, #number of localization windows
            W,# W chosen window length
            P,# P number of windows
            k,# k chosen modes
            a,
            e,
            t,
            data,            #T
            localized_data,  #loc_N x loc_len
            kept_RC,         # loc_N x loc_len x k
            EOF,     # W x k
            PC,      # k x P
            loc_emb,         # W x P
        )
    end


end

#----------------------------------------------------------------------
# helpers
#----------------------------------------------------------------------

#reconstruction ssa
function reconstructor(A,rho,N,M)
    P = N-M+1
    #A is the k_th PC (projection) of length P
    #rho is the k_th eigenvector (EOF) of length M

    #M is window length W
    R(t,M_t,L_t,U_t) = M_t*sum([A[t-j+1]*rho[j] for j in L_t:U_t])
    function choicer(t)
        if 1<=t<=M-1
            return 1/t,1,t
        elseif M<=t<=P
            return 1/M,1,M
        elseif P+1<=t<=N
            return 1/(N-t+1),t-N+M,M
        end
    end
    return [R(t,choicer(t)...) for t in 1:N]
end

#pricipal components
function pc!(X,rho,P,M)
    return unfold([Array(sum([X[t+j-1]*rho[j,:] for j in 1:M],dims = 1)...) for t in 1:P])
end
#make 2d from nested 2d array
unfold(multiArr) = reshape(hcat(multiArr...),size(multiArr)[1],size(multiArr[1])[1])

replicate(signal::Vector,O) = vcat(fill(signal,O)...)
retrieve_mid(signal::Vector,len) = signal[Int((length(signal)/2)-(len/2))+1:Int((length(signal)/2)+(len/2))]

#----------------------------------------------------------------------
# calculations
#----------------------------------------------------------------------
function c_emb(loc_window::Vector,N,W,P)
    Y = Float64[]
    for i=1:P
        Y = append!(Y,loc_window[i:i+W-1])
    end
    return Matrix(reshape(float.(Y),W,P)')
end
function c_EOF(emb_data::Matrix,k,a,e,t)     #returns the EOF[W,k]
    return Matrix(ManifoldLearning.transform(fit(DiffMap,emb_data,maxoutdim=k,t=t, α=a, ɛ=e))')
end
function c_PC(data::Vector,EOF::Matrix,P,W,k) #returns the PC[k,P],requires P that is size(emb)[1],requires W that is size(emb)[2]
    return Matrix(hcat([pc!(data,EOF[:,i],P,W) for i in 1:k]...)')
end
function c_RC(PC::Matrix,EOF::Matrix,loc_len,W,k)     #returns the RC[N,k]
    return Matrix(hcat([reconstructor(PC[i,:],EOF[:,i],loc_len,W) for i in 1:k]...))
end

#----------------------------------------------------------------------
# location window run
#----------------------------------------------------------------------
function loop_the_loc_w(self)

    loc_len = self.loc_len #len with repetitions
    loc_N = self.loc_N #number of localization windows
    W = self.W# W chosen window length
    P = self.P# P number of windows
    k = self.k# k chosen modes
    a = self.a
    e = self.e
    t = self.t

    for i = 1:loc_N
        self.loc_emb = c_emb(self.localized_data[i,:],loc_len,W,P)
        self.EOF = c_EOF(self.loc_emb,k,a,e,t)
        self.PC = c_PC(self.localized_data[i,:],self.EOF,P,W,k)
        self.kept_RC[i,:,:] = c_RC(self.PC,self.EOF,loc_len,W,k)
    end

    return self.kept_RC
end
#----------------------------------------------------------------------
# glue-together the mid-parts
#----------------------------------------------------------------------
glue(RC::Array) = hcat([retrieve_mid(RC) for i in 1:k]...)
#----------------------------------------------------------------------
# rescale and average
#----------------------------------------------------------------------
rescaler(loc_len) = [sin(pi*x/(loc_len-1)) for x=0:loc_len-1]

function averager(RC::Array,loc_N,loc_rate,W,k)
    list = push!([],zeros(W,k),[rescaler(loc_N) .* RC[i,:,:] for i in 1:loc_N]...,zeros(W,k))
    t=[]
    for i in 1:loc_N+1
        t = push!(t,list[i][1:Int(W*loc_rate),:].+list[i+1][Int(W*loc_rate)+1:end,:])
    end
return t
end
#----------------------------------------------------------------------
# parameter iteration?
#----------------------------------------------------------------------

#=
    


    #PREALLOCATION 

julia> function xinc!(ret::AbstractVector{T}, x::T) where T
           ret[1] = x
           ret[2] = x+1
           ret[3] = x+2
           nothing
       end;

julia> function loopinc_prealloc()
           ret = Vector{Int}(undef, 3)
           y = 0
           for i = 1:10^7
               xinc!(ret, i)
               y += ret[2]
           end
           return y
       end;

    function perform()
        signal::Vector,t,e,a,loc_size,N_repl,W,loc_rate,k=20

        P = Int(loc_size*N_repl - W + 1)
        locals_all = F1_locEmb(signal,loc_size,loc_rate,N_repl)
        println(size(locals_all))
        savepond = []
        for i in 1:size(locals_all)[1]
            locals = locals_all[i,:]
            emb = F2_emb(locals,W)
            EOF = F3_diff(emb,k,a,e,t)
            PC = F4_PC(locals,EOF,P,W,k)
            RC = F5_RC(PC,EOF,loc_size,W,N_repl,k)
            savepond = push!(savepond,RC) #EOF,PC
        end
        return vcat(savepond...) 
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

            =#
        