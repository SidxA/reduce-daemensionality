#load
function load_data()
    name = "testdata/FLX_US-PFa_FLUXNET2015_FULLSET_DD_1995-2014_1-4.csv"
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

    return df[:,20]
end

#centralize
function centralizer(data)
    m = mean(data)
    cv = std(data)
    return  (data.-m)./cv
end

#simple enough boundary hack for global ssa
function boundary_1(data)
    #first and last year are repeated
    #that means 365 datapoints
    first = data[1:365]
    last = data[end-364:end]
    return append!(first,data,last)
end

#simple enough cutoff
function cutoff(signal)
    return signal[366:end-365]
end

#embed lag
function embed_lag(data,W)
    Y = []
    for i=1:length(data)-W+1
        Y = append!(Y,data[i:i+W-1])
    end
    return reshape(float.(Y),W,length(data)-W+1)
end

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

#plot the ingredients
function plotthem(name,signal,EOFs,PC,RC,lambdas)
    #trajectory plot
    traj = scatter(signal,markersize=1,c=:blue,legend=:none,title="trajectory",alpha=0.5)
    traj =plot!(sum(RC,dims=2),c=:black)
    #spectrum
    spec = plot(lambdas,legend=:none,title="variance")
    #EOF
    eof = plot(EOFs,legend=:none,title="eof")
    #PC
    pc = plot(PC',legend=:none,title="pc")
    #RC
    rc = plot(RC,legend=:none,title="rc")
    savefig(plot(traj,spec,eof,pc,rc,layout=(5,1),dpi=800),dir*name)
    return nothing
end

#make 2d from nested 2d array
unfold(multiArr) = reshape(hcat(multiArr...),size(multiArr)[1],size(multiArr[1])[1])

# windows at W/2 along the data cuts off the rest
function embed_lag_centered(x,win_l)
    T = length(x)
    P = win_l
    N = floor(Int,(2(T+1)/P))-1
    data = zeros(N,P)
    ind = 1
    for t in 1:N
        ind_l = ind + P - 1
        if ind_l <= T
            data[t,:] = x[ind:ind_l]
            ind += Int(P/2)
        else
            break
        end
    end
    return data
end

#simple enough boundary hack for local ssa
boundary_replicate(data) = append!([vcat(data) for i=1:5]...)

#----------------------------------------------------------------------
# global SSA
#----------------------------------------------------------------------
#returns rc, possiblysignal,eof,pc,rc,lambdas
function SSA_g(data,W,key=:auto)
    k=20 #hardcoded kept number of EOF
    data = data |> boundary_1 |> centralizer

    N=length(data)
    P = N - W +1

    emb_data = embed_lag(data,W)
    pca = fit(PCA,emb_data,method=key,pratio=1.0,maxoutdim=k)

    EOFs = projection(pca)
    PC = hcat([pc!(data,EOFs[:,i],P,W) for i in 1:k]...)'
    RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) |> cutoff for i in 1:k]...)

    return cutoff(data),EOFs,PC,RC,principalvars(pca)
end

#----------------------------------------------------------------------
# local SSA
#----------------------------------------------------------------------
#returns rc (possibly signal)
function local_SSA(signal,W,P)
    k=20 #hardcoded number of modes

    windowed_data = embed_lag_centered(signal,W) # size W x floornumber
    rescaler = [sin(pi*x/(W-1)) for x=0:W-1]

    floornumber = size(windowed_data)[1]

    local_rcs = []
    pca = PCA
    data = zeros(5*W)

    #given one window
    #replicate it five times
    #perform the ssa - embed (1.25 a)
    #now take the EOF of this
    #we rescale them (sin of length P/2)

    for window_number in 1:floornumber
        data = windowed_data[window_number,:] |> boundary_replicate |> centralizer
        N = length(data)
        K = N - P + 1

        emb_data = embed_lag(data,P)
        pca = fit(PCA,emb_data,maxoutdim=k,pratio=1.0,method=:svd)
        EOFs = projection(pca)

        PC = hcat([pc!(data,EOFs[:,i],K,P) for i in 1:k]...)'

        RC = hcat([reconstructor(PC[i,:],EOFs[:,i],W,P) for i in 1:k]...)

        local_rcs = push!(local_rcs,rescaler.*RC) #[Int(4*W/2):Int(6*W/2)-1,:]
    end
    #we add them together at the mid-endpoints
    #list together with zero window on beginning and end

        list= []
        list = push!(list,zeros(W,k),local_rcs...,zeros(W,k))

        t=[]
        for i in 1:floornumber+1
            t = push!(t,list[i][1:Int(W/2),:].+list[i+1][Int(W/2)+1:end,:])
        end

    #now we have a list of floornumber x wind/2 x k
    rc_glob = vcat([t[i] for i=1:length(t)]...)
    return rc_glob#signal[1:size(rc_glob)[1]],rc_glob
end

#----------------------------------------------------------------------
# global Diffusion map
#----------------------------------------------------------------------
#returns rc
function global_diffmap(signal,W)
    k=20 #hardcoded kept number of EOF
    data = signal |> boundary_1 |> centralizer

    N=length(data)
    P = N - W +1

    diff = fit(DiffMap,embed_lag(data,W)',maxoutdim=k,t=1, α=0.5, ɛ=512)

    EOFs = ManifoldLearning.transform(diff)'
    PC = hcat([pc!(data,EOFs[:,i],P,W) for i in 1:k]...)'
    RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) |> cutoff for i in 1:k]...)

    return RC
end

#----------------------------------------------------------------------
# local Diffusion map
#----------------------------------------------------------------------
#returns RC... i guess
replicate(signal::Vector,O) = vcat(fill(signal,O)...)
retrieve_mid(signal::Vector,len) = signal[Int((length(signal)/2)-(len/2))+1:Int((length(signal)/2)+(len/2))]
#simple replication of the local window a number of times, put odd number of replication!
function F1_locEmb(signal::Vector,loc_size,loc_rate,N_repl)
    @assert isodd(N_repl) "replication not odd"
    signal = isodd(length(signal)) ? signal[1:end-1] : signal
    T = length(signal)
    N_locs = floor(Int,((T+1)/(loc_size * loc_rate)))
    data = zeros(N_locs,loc_size*N_repl)
    println("replicate loc_emb of signal_l $T")
    println("windows $loc_size, cover_rate $loc_rate, N_repl $N_repl, throwing away $(Int(T-N_locs*loc_size*loc_rate))")
    ind = 1
    for t in 1:N_locs
        ind_l = ind + loc_size - 1 #position_end is position + loc_size -1
        if ind_l < T #if position_end in the series
            data[t,:] = replicate(signal[ind:ind_l],N_repl) #put the window
            ind += Int(loc_size * loc_rate) #position shifts by loc_size / loc_rate
        else
            break
        end
    end
    return data
end

#this is not centralizing the data, maybe do for comparison!
function F2_emb(loc_window::Vector,W)
    N=length(loc_window)
    @assert ((N+1)/2>W) "choose smaller W"
    P = N - W +1
    Y = []
    for i=1:length(loc_window)-W+1
        Y = append!(Y,loc_window[i:i+W-1])
    end
    return Matrix(reshape(float.(Y),W,length(loc_window)-W+1)')
end

#returns the EOF[W,k]
function F3_diff(emb_data::Matrix,k=20,a=0.5,e=512,t=1)
        return Matrix(ManifoldLearning.transform(fit(DiffMap,emb_data,maxoutdim=k,t=t, α=a, ɛ=e))')
end

#returns the PC[k,P],requires number of embedding windows P that is size(emb)[1],requires embedding window size that is size(emb)[2]
function F4_PC(data::Vector,EOF::Matrix,P,W,k=20)
    return Matrix(hcat([pc!(data,EOF[:,i],P,W) for i in 1:k]...)')
end

# this we could also compare to something thtat does nto cut off but rather averages
#returns the RC[N,k]
function F5_RC(PC::Matrix,EOF::Matrix,loc_size,W,N_repl,k=20)
    return hcat([retrieve_mid(reconstructor(PC[i,:],EOF[:,i],loc_size*N_repl,W),loc_size) for i in 1:k]...)
end

function l_diff(signal::Vector,loc_size,loc_rate,N_repl,W,k=20)

    P = Int(loc_size*N_repl - W + 1)
    locals_all = F1_locEmb(signal,loc_size,loc_rate,N_repl)
    println(size(locals_all))
    savepond = []
    for i in 1:size(locals_all)[1]
        locals = locals_all[i,:]
        emb = F2_emb(locals,W)
        EOF = F3_diff(emb,k)
        PC = F4_PC(locals,EOF,P,W,k)
        RC = F5_RC(PC,EOF,loc_size,W,N_repl,k)
        savepond = push!(savepond,RC) #EOF,PC
    end
    return vcat(savepond...) 
end

#----------------------------------------------------------------------
# simple hilbert transform via FourierAnalysis.analyticsignal
#----------------------------------------------------------------------

#phase transformation
phases(an_s::Vector) = atan.(imag.(an_s)./real.(an_s))

function iterate_hilb(an_s::Vector,k)
    for i=1:k
        an_s = phases(analyticsignal(an_s))
    end
    return an_s
end

#count sign flips
count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])

#count maxima
count_maxima(signal) = length(findmaxima(signal)[1])




#estimating of period length via phase estimation#max should be good enough
#only hilbert trafo
function phase_information(phases) #faulty!
    ind,val = findmaxima(phases)
    T = length(phases)/length(ind)
    outliers = length(filter(x -> (abs(x-pi/2)>0.1), val))/length(ind)
    return T,val,outliers
end


#protophase
function protophase(signal::Vector)
    x = signal
    y = imag.(hilbert_transform(signal))
    curvelength = Array{Float64}(undef,0)
    curvelength = push!(curvelength,0.0)
    for index = 2:length(signal)
        distance = sqrt((x[index]-x[index-1])^2 + (y[index]-y[index-1])^2)
        curvelength = push!(curvelength,curvelength[index-1]+distance)
    end
    return x,y,curvelength
end

#cutoff for the intervals for a protophase going from 0...N to range alwayss between 0 and 2pi and again...
function crude_normalization(protophase::Vector)
    for (i,x) in enumerate(protophase)
        if x >= 2pi
            protophase[i:end] .-= 2*pi
        end
    end
    return protophase
end

#----------------------------------------------------------------------
# local SSA
#----------------------------------------------------------------------
#Not working?

#returns rc (possibly signal)
function local_SSA(signal,W,P)
    k=20 #hardcoded number of modes

    windowed_data = embed_lag_centered(signal,W) # size W x floornumber
    rescaler = [sin(pi*x/(W-1)) for x=0:W-1]

    floornumber = size(windowed_data)[1]

    local_rcs = []
    pca = PCA
    data = zeros(5*W)

    #given one window
    #replicate it five times
    #perform the ssa - embed (1.25 a)
    #now take the EOF of this
    #we rescale them (sin of length P/2)

    for window_number in 1:floornumber
        data = windowed_data[window_number,:] |> boundary_replicate |> centralizer
        N = length(data)
        K = N - P + 1

        emb_data = embed_lag(data,P)
        pca = fit(PCA,emb_data,maxoutdim=k,pratio=1.0,method=:svd)
        EOFs = projection(pca)

        PC = hcat([pc!(data,EOFs[:,i],K,P) for i in 1:k]...)'

        RC = hcat([reconstructor(PC[i,:],EOFs[:,i],W,P) for i in 1:k]...)

        local_rcs = push!(local_rcs,rescaler.*RC) #[Int(4*W/2):Int(6*W/2)-1,:]
    end
    #we add them together at the mid-endpoints
    #list together with zero window on beginning and end

        list= []
        list = push!(list,zeros(W,k),local_rcs...,zeros(W,k))

        t=[]
        for i in 1:floornumber+1
            t = push!(t,list[i][1:Int(W/2),:].+list[i+1][Int(W/2)+1:end,:])
        end

    #now we have a list of floornumber x wind/2 x k
    rc_glob = vcat([t[i] for i=1:length(t)]...)
    return rc_glob#signal[1:size(rc_glob)[1]],rc_glob
end


#----------------------------------------------------------------------
# local Diffusion map
#----------------------------------------------------------------------
#working but needs sorting out for the replications?
#returns RC... i guess
replicate(signal::Vector,O) = vcat(fill(signal,O)...)
retrieve_mid(signal::Vector,len) = signal[Int((length(signal)/2)-(len/2))+1:Int((length(signal)/2)+(len/2))]
#simple replication of the local window a number of times, put odd number of replication!
function F1_locEmb(signal::Vector,loc_size,loc_rate,N_repl)
    @assert isodd(N_repl) "replication not odd"
    signal = isodd(length(signal)) ? signal[1:end-1] : signal
    T = length(signal)
    N_locs = floor(Int,((T+1)/(loc_size * loc_rate)))
    data = zeros(N_locs,loc_size*N_repl)
    #println("replicate loc_emb of signal_l $T")
    #println("windows $loc_size, cover_rate $loc_rate, N_repl $N_repl, throwing away $(Int(T-N_locs*loc_size*loc_rate))")
    ind = 1
    for t in 1:N_locs
        ind_l = ind + loc_size - 1 #position_end is position + loc_size -1
        if ind_l < T #if position_end in the series
            data[t,:] = replicate(signal[ind:ind_l],N_repl) #put the window
            ind += Int(loc_size * loc_rate) #position shifts by loc_size / loc_rate
        else
            break
        end
    end
    return data
end

#this is not centralizing the data, maybe do for comparison!
function F2_emb(loc_window::Vector,W)
    N=length(loc_window)
    @assert ((N+1)/2>W) "choose smaller W"
    P = N - W +1
    Y = []
    for i=1:length(loc_window)-W+1
        Y = append!(Y,loc_window[i:i+W-1])
    end
    return Matrix(reshape(float.(Y),W,length(loc_window)-W+1)')
end

#returns the EOF[W,k]
function F3_diff(emb_data::Matrix,k=20,a=0.5,e=512,t=1)
        return Matrix(ManifoldLearning.transform(fit(DiffMap,emb_data,maxoutdim=k,t=t, α=a, ɛ=e))')
end

#returns the PC[k,P],requires number of embedding windows P that is size(emb)[1],requires embedding window size that is size(emb)[2]
function F4_PC(data::Vector,EOF::Matrix,P,W,k=20)
    return Matrix(hcat([pc!(data,EOF[:,i],P,W) for i in 1:k]...)')
end

# this we could also compare to something thtat does nto cut off but rather averages
#returns the RC[N,k]
function F5_RC(PC::Matrix,EOF::Matrix,loc_size,W,N_repl,k=20)
    return hcat([retrieve_mid(reconstructor(PC[i,:],EOF[:,i],loc_size*N_repl,W),loc_size) for i in 1:k]...)
end

function l_diff(signal::Vector,t,e,a,loc_size,N_repl,W,loc_rate,k=20)

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


#plot the ingredients
function plotthem(name,signal,EOFs,PC,RC,lambdas)
    #trajectory plot
    traj = scatter(signal,markersize=1,c=:blue,legend=:none,title="trajectory",alpha=0.5)
    traj =plot!(sum(RC,dims=2),c=:black)
    #spectrum
    spec = plot(lambdas,legend=:none,title="variance")
    #EOF
    eof = plot(EOFs,legend=:none,title="eof")
    #PC
    pc = plot(PC',legend=:none,title="pc")
    #RC
    rc = plot(RC,legend=:none,title="rc")
    savefig(plot(traj,spec,eof,pc,rc,layout=(5,1),dpi=800),dir*name)
    return nothing
end

function embed_lag(x,win_l)
    T = length(x)
    P = win_l
    N = floor(Int,(2(T+1)/P))-1
    data = zeros(N,P)
    ind = 1
    for t in 1:N
        ind_l = ind + P - 1
        if ind_l < T
            data[t,:] = x[ind:ind_l]
            ind += Int(P/2)
        else
            break
        end
    end
    return data
end

function dothepca(c_data,window_length,dims)
    X = convert(Matrix, embed_lag(c_data,window_length))
    pca = fit(PCA,X,method=:auto,maxoutdim=dims,pratio=1.0)
    t = transform(pca,X)
    z = zeros(dims)
    p = principalvars(pca)
    l = length(p)
    z[1:l] = p 
    return z 
end

function loopthepca(c_data,length_array,dims)
    l = length(length_array)
    z = zeros(l,dims)
    for (i,x) in enumerate(length_array)
        z[i,:] = dothepca(c_data,x,dims)
    end
    return z 
end

function centralize(data)
    m = mean(data)
    cv = std(data)
    return  (data.-m)./cv
end

function deseasonalize()
    for name in nameslist
        signal = npzread("/net/scratch/lschulz/cropdata/"*name)
        N = length(signal)
        P = N-W+1
        pca = fit(PCA,embed_lag(signal,W),method=:svd,pratio=1.0,maxoutdim=k)
        EOF = projection(pca)
        PC = hcat([pc!(signal,EOF[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOF[:,i],N,W) for i in 1:k]...)
        npzwrite("/net/scratch/lschulz/cropdata_deseason/"*name, signal .- sum(RC[:,1:2],dims=2))
        println(name)
    end
end