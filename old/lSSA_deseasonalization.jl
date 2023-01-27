#repeated windows at W/2 along the data cuts off the rest
function repetition_windows(x::Vector{Float32},Rep_len,Rep_num)
    N = length(x)

    floornumber = floor(Int,(2(N+1)/Rep_len))-1

    data = zeros(Rep_len*Rep_num,floornumber)

    ind = 1

    for t in 1:floornumber
        ind_l = ind + Rep_len - 1
        if ind_l <= N
            data[:,t] = repeat(x[ind:ind_l],Rep_num)
            ind += Int(Rep_len/2)
        else
            break
        end
    end
    return data
end

function centralized_embed_lag(data::Vector{Float32},W::Int64)
    #centralize
    function centralizer(data)
        m = mean(data)
        cv = std(data)
        return  (data.-m)./cv
    end
    Y = []

    for i=1:length(data)-W+1

        Y = append!(Y,centralizer(data[i:i+W-1]))

    end
    return reshape(float.(Y),W,length(data)-W+1)::Matrix{Float32}
end

retrieve_mid(signal::Vector,len) = signal[Int((length(signal)/2)-(len/2))+1:Int((length(signal)/2)+(len/2))]

#----------------------------------------------------------------------
# local SSA
#----------------------------------------------------------------------
#only intended for deseasonalization, build simple for fixed 5* 1 year with 1.25 embedding
mutable struct deseason_lSSA
    N::Int64
    k::Int64
    signal::Vector{Float32}

    #window size 1 year
    Rep_len::Int64
    #replication number 5
    Rep_num::Int64
    #embedding inside 1,2 years
    Rep_W::Int64
    #number of repetition windows
    Rep_W_num::Int64

    Rep_M::Matrix{Float32}
    Rep_emb_M::Matrix{Float32}
    
    EOF::Vector{Float32}
    PC::Vector{Float32}
    RC::Vector{Float32}

    RC_save::Matrix{Float32}

    function deseason_lSSA(signal::Vector{Float32},Rep_W)
        N               = length(signal)
        k               = 2
        Rep_len         = 366 #needs to be an even number (midpoints...)
        Rep_num         = 5
        #Rep_W           = 438 #why not do this for 365
        Rep_W_num       = floor(Int,(2(N+1)/Rep_len))-1

        Rep_M           = repetition_windows(signal,Rep_len,Rep_num)
        Rep_emb_M       = Array{Float32}(undef,(Rep_len*Rep_num)-Rep_W+1,Rep_W)

        #cutoff: give new signal and new N
        N               = Rep_len*Int(floor(N/Rep_len))
        signal          = signal[1:N]
        EOF             = Vector{Float32}(undef,Rep_W)
        PC              = Vector{Float32}(undef,(Rep_len*Rep_num)-Rep_W+1)
        RC              = Vector{Float32}(undef,(Rep_len*Rep_num))

        RC_save         = Array{Float32}(undef,Rep_len,Rep_W_num)
        return new(
            N,
            k,
            signal,
            Rep_len,
            Rep_num,
            Rep_W,
            Rep_W_num,
            Rep_M,
            Rep_emb_M,
            EOF,
            PC,
            RC,
            RC_save
        )
    end
end
    #fixed N, W and k, list of loc
    #then iterates the parameters
    #with centralization!
    #saves as name number N W k 

function deseasonalize_lSSA(signal::Vector{Float32},sampleyear)

    l = deseason_lSSA(signal,sampleyear)

    N = l.Rep_len*l.Rep_num
    W = l.Rep_W
    P = N - W +1
    rescaler = [sin(pi*x/(W-1)) for x=0:l.Rep_len-1]
    #given one window
    #replicate it five times
    #perform the ssa - embed (1.2 a)
    #now take the EOF of this
    #we rescale them (sin of length P/2)
    for win_i in 1:l.Rep_W_num
        l.Rep_emb_M = centralized_embed_lag(l.Rep_M[:,win_i],l.Rep_W)
        l.EOF = sum(projection(fit(PCA,l.Rep_emb_M::Matrix{Float32},maxoutdim=l.k)),dims=2)[:,1]
        l.PC = pc!(l.Rep_M[:,win_i],l.EOF,P,W)[:,1]
        l.RC = Float32.(reconstructor(l.PC,l.EOF,N,W))
        l.RC_save[:,win_i] = retrieve_mid(l.RC,l.Rep_len) .* rescaler
    end
    #list together with zero window on beginning and end
    list = hcat(zeros(l.Rep_len),[l.RC_save[:,i] for i in 1:l.Rep_W_num]...,zeros(l.Rep_len))
    season=Float32[]
    #we add them together at the mid-endpoints
    for i in 1:l.Rep_W_num+1
        season = append!(season,Float32.(list[Int(l.Rep_len/2)+1:end,i].+list[1:Int(l.Rep_len/2),i+1]))
    end
    #return signal, signal .- season, var(season) / var(signal)
    return Float32.(signal.-vcat(season,rand((-1/2:0.01:1/2),length(signal)-length(season))))
end

function deseasonalize_gSSA(signal::Vector{Float32},sampleyear)
    N = length(signal)
    W = sampleyear
    P = N-W+1
    k=2

    emb = centralized_embed_lag(signal,W)
    EOF = sum(projection(fit(PCA,emb,method=:svd,pratio=1.0,maxoutdim=k)),dims=2)
    PC = (EOF' * emb)[1,:]
    #PC = pc!(signal,EOF,P,W)[:,1]
    season = reconstructor(PC,EOF,N,W)
    return Float32.(signal .- season)
end