using Distributed
@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")

    function sample_par!()

        #random data set from the flux sites
        Lnames = readdir("/net/home/lschulz/logs/KW_26/TA_ERA_DAY")
        name = rand(Lnames)
        #random boundary extension
        B_L = Int64[0,Int.(floor.(91.3:91.3:732).+1)...]
        B = rand(B_L)
        # from dataset with boundary, know max W and then choose random W
        data = boundary(centralizer(npzread("/net/home/lschulz/logs/KW_26/TA_ERA_DAY/"*name)),B)
        N = length(data)
        W_L = Int64[Int.(floor.(730.5:182.625:N/2-2).+1)...]
        W = rand(W_L)
        #random number of k vectors
        K_L = Int64[16:4:64...]
        k=rand(K_L)
        #random diffusion parameters
        alpha_L = Float64[0.5,1.0]
        alpha = rand(alpha_L)
        t_L = Int64[1:1:5...]
        t = rand(t_L)
        eps_L = Int64[12:12:200...]
        eps = rand(eps_L)

        return name,data,B,N,W,k,alpha,t,eps
    end

    function sampler!()
        name,signal,B,N,W,k,alpha,t,eps = sample_par!()
        P=N-W+1
        savedir = dir*"sample/$(name[1:6])_B$(B)_W$(W)_k$(k)_a$(alpha)_t$(t)_e$(eps)"

        #EOF_SSA = projection(fit(PCA,embed_lag(signal,W),method=:auto,pratio=1.0,maxoutdim=k))
        EOF = ManifoldLearning.transform(fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=alpha, ɛ=eps))'

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
        npzwrite(savedir*"_Diff",EOF)

        return nothing #result
    end



    #logfile = open(dir*"samplelog.txt", "a") 
end


Threads.@threads for i=1:10000
    sampler!()
    #write(logfile, "sample \n")
    println("sample")
    if time()>1.65701511097627e9 break end
end

#@everywhere close(logfile)

