#first lets do some fourier
function fourier(data)
    signal = boundary_1(load_data())
    # Number of points 
    N = length(signal)
    # Sample period
    Ts = 1/365
    # Start time 
    t0 = 0
    tmax = t0 + (N-1) * Ts
    # time coordinate
    t = t0:Ts:tmax


    # Fourier Transform of it 
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    #years = [i/365 for i in 1:8]
    # plots 
    savefig(plot(t, signal, title = "Signal"),dir*"signal")
    savefig(plot(freqs, abs.(F), title = "Spectrum",ylim=(0,10000),xlim=(0,200)),dir*"spectrum")

end
#lets do some global ssa returns eof,pc,rc
function SSA_g(data,W)
    #hardcod
    k=20

    data = data |> boundary_1 |> centralizer
    N=length(data)
    P = N - W +1

    #embedding matrix / maybe already do the centered one?
    emb_data = embed_lag(data,W)
    pca = fit(PCA,emb_data,maxoutdim=k,pratio=1.0,method=:svd)

    EOFs = projection(pca)
    PC = hcat([pc!(data,EOFs[:,i],P,W) for i in 1:k]...)'
    RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) |> cutoff for i in 1:k]...)

    return EOFs,PC,RC
end
#reconstruction
function transform_global_SSA(EOF,PC,N,W,P)

    #EOF = vec(projection(pca)[:,season])#we have the EOF of length W
    #PC = vec(MultivariateStats.transform(pca,data)[season,:])#we have the PC of length P=N-W+1
    #we have the lagged trajectory of size (W,P)

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
        RC = push!(RC,(1/M)*sum(inner_sum))
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
        RC = push!(RC,(1/M)*sum(inner_sum))
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
        RC = push!(RC,(1/M)*sum(inner_sum))
    end

    return RC
end

function this()
    #parameter test for th ediffusion map => epsilon was good around 2-500 !
    for (i,epsilon)=enumerate(2 .^Array(-6.0:1.0:8.0)), W in 200:200:1000
        savefig(plot(ManifoldLearning.transform(fit(DiffMap,embed_lag(signal,W),
        maxoutdim=10,t=1, α=1.0, ɛ=epsilon))',
        legend=:none,title=string(epsilon," ",W)),
        dirdiff*string("_",W,"_",i))
    end
end


function the_whole_run(signal,W)
    data = centralizer(signal)
    N_hilbert = 20
    hilbert_data = data
    for i in 1:N
        hilbert_data = phases(analyticsignal(hilbert_data))
    end
    data = embed_lag(data,W)
    components_P = ManifoldLearning.transform(fit(DiffMap,data,maxoutdim=100,t=1, α=0.5, ɛ=512))

end


#=
#for W=Array(400:400:2000), boundary_cutoff = Array(400:400:2000)
for W=Array(200:200:4000), boundary_cutoff = Array(400:400:4000)
    #simple enough boundary hack for global ssa
    boundary_1(data) =  append!(data[1:boundary_cutoff],data,data[end-boundary_cutoff:end])

    #simple enough cutoff
    cutoff(signal) =signal[boundary_cutoff+1:end-boundary_cutoff]

    npzwrite(dir*"runs/"*string("RC_components_W_",W,"_B_",boundary_cutoff),global_diffmap(load_data(),W))
end

W = 3200
B = 4000



=#

function SSA_g2(data,W,key=:auto)
    k=20 #hardcoded kept number of EOF
    data = data |> boundary_1 |> centralizer

    N=length(data)
    P = N - W +1

    emb_data = embed_lag(data,W)
    pca = fit(PCA,emb_data,method=key,pratio=1.0,maxoutdim=k)

    EOFs = projection(pca)
    PC = hcat([pc!(data,EOFs[:,i],P,W)  for i in 1:k]...)'
    RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) |> cutoff for i in 1:k]...)

    return pca
end


#=
for W=Array(200:200:4000), boundary_cutoff = Array(400:400:4000)
    #simple enough boundary hack for global ssa
    boundary_1(data) =  append!(data[1:boundary_cutoff],data,data[end-boundary_cutoff:end])

    #simple enough cutoff
    cutoff(signal) =signal[boundary_cutoff+1:end-boundary_cutoff]

    npzwrite(dir*"runs/"*string("RC_SSA_g_W_",W,"_B_",boundary_cutoff),SSA_g2(load_data(),W))
end

=#
function F1_locEmb_plain(signal::Vector,loc_size,loc_rate)
    T = length(signal)
    N_locs = floor(Int,((T+1)/(loc_size * loc_rate)))
    data = zeros(N_locs,loc_size)
    println("plain loc_emb of signal_l $T")
    println("windows $loc_size, cover_rate $loc_rate, throwing away $(Int(T-N_locs*loc_size*loc_rate))")

    ind = 1
    for t in 1:N_locs
        ind_l = ind + loc_size - 1 #position_end is position + loc_size -1
        if ind_l < T #if position_end in the series
            data[t,:] = signal[ind:ind_l] #put the window
            ind += Int(loc_size * loc_rate) #position shifts by loc_size / loc_rate
        else
            break
        end
    end
    return data
end




#EOFs = ManifoldLearning.transform(diff)'
#PC = hcat([pc!(data,EOFs[:,i],P,W) for i in 1:k]...)'
#RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) |> cutoff for i in 1:k]...)


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
function F5_RC_cutoff(PC::Matrix,EOF::Matrix,loc_size,W,N_repl,k=20)
    return hcat([retrieve_mid(reconstructor(PC[i,:],EOF[:,i],loc_size*N_repl,W),loc_size) for i in 1:k]...)
end

rescaler(N_loc) = [sin(pi*x/(N_loc-1)) for x=0:N_loc-1]

#we need to do this via avergaing
function F5_RC_averaging(PC::Matrix,EOF::Matrix,loc_size,W,N_repl,k=20)

    rescaler.*reconstructor(PC[i,:],EOF[:,i],loc_size*N_repl,W) # for j in N_locs, i in k
    
    list= []
    list = push!(list,zeros(W,k),local_rcs...,zeros(W,k))

    t=[]
    for i in 1:floornumber+1
        t = push!(t,list[i][1:Int(W/2),:].+list[i+1][Int(W/2)+1:end,:])
    end

#now we have a list of floornumber x wind/2 x k
rc_glob = vcat([t[i] for i=1:length(t)]...)


    return hcat([
        retrieve_mid(
            reconstructor(PC[i,:],EOF[:,i],loc_size*N_repl,W)
            ,loc_size)
             for i in 1:k]
             ...)
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
