using Distributed
@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")

    #return name,data,B,N,Ws,ks,W,k
    function sample_ssa_par!()

        #random data set from the flux sites
        Lnames = readdir("/net/home/lschulz/logs/KW_26/TA_ERA_DAY")
        name = rand(Lnames)
        #random boundary extension
            #B_L = Int64[0,Int.(floor.(91.3:91.3:732).+1)...]
        B = rand([182,365,365+182])#rand(B_L)
        # from dataset with boundary, know max W and then choose random W
        data = boundary(centralizer(npzread("/net/home/lschulz/logs/KW_26/TA_ERA_DAY/"*name)),B)
        N = length(data)
        #W for deseasoning
            #Ws_L = Int64[Int.(floor.(182.625:182.625:730.5).+1)...]
        Ws = rand([365,547])#rand(Ws_L)
        #random number of k vectors for deseasoning
            #Ks_L = Int64[4:4:16...]
        ks=12#rand(Ks_L)
        #long W for modes
        W_L = Int64[Int.(floor.((365.25*7):182.625:N/2).+1)...]
        W = rand(W_L)
        #random number of k vectors
            #K_L = Int64[16:4:64...]
        k=32#rand(K_L)

        return name,data,B,N,Ws,ks,W,k
    end

    #return signal,season
    function de_seasonalizer!(signal,Ws,ks)
        W=Ws
        k=ks
        N = length(signal)
        P = N-W+1
        pca = fit(PCA,embed_lag(signal,W),method=:svd,pratio=1.0,maxoutdim=k)
        EOF = projection(pca)
        PC = hcat([pc!(signal,EOF[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOF[:,i],N,W) for i in 1:k]...)
        return sum(RC[:,1:2],dims=2)
    end


    function sampler_ssa!()
        name,data,B,N,Ws,ks,W,k = sample_ssa_par!()
        title = "$(name[1:6])_B$(B)_W$(W)_k$(k)_Ws$(Ws)_ks$(ks)"
        savedir = dir*"sample_ssa2/"*title

        season = de_seasonalizer!(data,Ws,ks)
        signal = data .- season
        N = length(signal)
        P = N-W+1
        pca = fit(PCA,embed_lag(signal,W),method=:svd,pratio=1.0,maxoutdim=k)
        EOF = projection(pca)
        npzwrite(savedir,EOF)
        npzwrite(savedir*"_signal",hcat(data,signal,season))
        println(title)
        #PC = hcat([pc!(signal,EOF[:,i],P,W) for i in 1:k]...)'
        #RC = hcat([reconstructor(PC[i,:],EOF[:,i],N,W) for i in 1:k]...)
        #EOF_SSA = projection(fit(PCA,embed_lag(signal,W),method=:auto,pratio=1.0,maxoutdim=k))
        #EOF = ManifoldLearning.transform(fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=alpha, ɛ=eps))'

        #epsi_T = 30
        #epsi_A = 0.1
        #per_eps = 350

        #pairs,periods = pairing_full(EOF,epsi_T,epsi_A,false)
        #per_cc = periods[periods .> per_eps]
        #pe = Float64[0.5 .* N ./ iterated_hilbert(EOF[:,i],64)[1] for i=1:k]
        #per_hht = pe[pe .> per_eps]
        #result = ""
        #result *= "$(name[1:6]) \t $(B) \t $(W) \t $(k) \t $(alpha) \t $(t) \t $(eps) \n"
        #result *= "Diff CC $(per_cc) \t HHT $(per_hht) \n"

        return nothing #result
    end

end


Threads.@threads for i=1:100000000
    sampler_ssa!()
    #write(logfile, "sample \n")
    if time()>1.657083507354178e9 break end
end

include("/net/home/lschulz/scripts/sampling_ssa_read.jl")

