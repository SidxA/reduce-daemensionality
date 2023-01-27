using Distributed
@everywhere begin
include("/net/home/lschulz/scripts/toolbox.jl")

#----------------------------------------------------------------------
# local SSA for deseason WITH CHoice of Rep-W
#----------------------------------------------------------------------
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

    function deseason_lSSA(signal::Vector{Float32},Rep_W::Int64)
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

#return season
function deseasonalize_lSSA(signal::Vector{Float32},Rep_W::Int64)

    l = deseason_lSSA(signal,Rep_W)

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

    return season
end

#----------------------------------------------------------------------
# global SSA for deseason WITH W = 366
#----------------------------------------------------------------------
#return season
function deseasonalize_gSSA(signal::Vector{Float32})
    N = length(signal)
    W = 366
    P = N-W+1
    k=2

    emb = centralized_embed_lag(signal,W)
    EOF = sum(projection(fit(PCA,emb,method=:svd,pratio=1.0,maxoutdim=k)),dims=2)[:,1]
    PC = pc!(signal,EOF,P,W)[:,1]
    season = reconstructor(PC,EOF,N,W)
    return season
end
#----------------------------------------------------------------------
# global SSA 
#----------------------------------------------------------------------
mutable struct matrix_holder
    #the variables
    N::Int64
    W::Int64
    k::Int64
    signal::Vector{Float32}
    emb::Matrix{Float32}
    EOF::Matrix{Float32}
    PC::Matrix{Float32}
    RC::Matrix{Float32}
    #the init function
    function matrix_holder(N::Int,W::Int,k::Int,signal::Vector{Float32})
        emb             = Array{Float32}(undef,N-W+1,W) #P x W
        EOF             = Array{Float32}(undef,W,k)
        PC             = Array{Float32}(undef,k,N-W+1)
        RC             = Array{Float32}(undef,N,k)
        return new(N,W,k,signal,emb,EOF,PC,RC)
    end
end

#----------------------------------------------------------------------
# one package calculation takes care of one location
#Boundary yes no
#W is the max W is the max fullyear
#with 4 lSSA and 1 gSSA this gives rise to 5*2*2 = 20 calculations per package
#----------------------------------------------------------------------
function SSA_package_handler(data::Matrix{Float32},package,outdir::String,nameslist::Vector{String})

    loc_ind = package[1]::Int64
    data_ind = package[2]::Int64
    name = nameslist[loc_ind]
    k = 32

    for B = [true,false], W_f = [true,false]
        signal = Float32.(B ? boundary(data[:,data_ind],365) : data[:,data_ind])

        for Rep_W in [456,438,366,416]
            season = Float32.(deseasonalize_lSSA(signal,Rep_W))
            signal = signal[1:length(season)] .- season
            N = length(signal)

            W = Int(iseven(Int(floor(N/2 - 1 ))) ? Int(floor(N/2 - 1 )) : Int(floor(N/2 - 2 )))
            W = W_f ? Int(floor(W-W%365.25)) : W

            P = N - W +1

            object = matrix_holder(N,W,k,signal::Vector{Float32})
            object.emb = centralized_embed_lag(object.signal,W)'
            object.EOF = Matrix(projection(fit(PCA,object.emb',method=:svd,pratio=1.0,maxoutdim=k)))
            object.PC = hcat([pc!(object.signal,object.EOF[:,i],P,W) for i in 1:k]...)'
            object.RC = hcat([reconstructor(object.PC[i,:],object.EOF[:,i],N,W) for i in 1:k]...)

            B ? object.RC = hcat([cutoff(object.RC[i,:],365) for i=1:k]...)  : object.RC = object.RC

            jldsave(outdir*"SSA_$(loc_ind)_$(name)_$(N)_$(W)_$(k)_$(B)_$(W_f)_lSSA_$(Rep_W).jld2",
            EOF = Matrix(object.EOF),PC = Matrix(object.PC),RC = Matrix(object.RC))
        end

        season = Float32.(deseasonalize_gSSA(signal))
        signal = signal[1:length(season)] .- season
        N = length(signal)

        W = Int(iseven(Int(floor(N/2 - 1 ))) ? Int(floor(N/2 - 1 )) : Int(floor(N/2 - 2 )))
        W = W_f ? Int(floor(W-W%365.25)) : W

        P = N - W +1

        object = matrix_holder(N,W,k,signal)
        object.emb = centralized_embed_lag(object.signal,W)'
        object.EOF = Matrix(projection(fit(PCA,object.emb',method=:svd,pratio=1.0,maxoutdim=k)))
        object.PC = hcat([pc!(object.signal,object.EOF[:,i],P,W) for i in 1:k]...)'
        object.RC = hcat([reconstructor(object.PC[i,:],object.EOF[:,i],N,W) for i in 1:k]...)

        object.RC = Float32.(B ? hcat([cutoff(object.RC[i,:],365) for i=1:k]...) : object.RC)

        jldsave(outdir*"SSA_$(loc_ind)_$(name)_$(N)_$(W)_$(k)_$(W_f)_$(B)_gSSA_0.jld2",
        EOF = Matrix(object.EOF),PC = Matrix(object.PC),RC = Matrix(object.RC))
    end
    return nothing
end
# for packages we basically only need the choice of the location
#right now the packs only has 1 entry with loc_ind
function packs(spotchoice)
    p = []
    for (i,x) in enumerate(spotchoice)
        p = push!(p,(x,i))
    end
    return p
end



#need to fill in the new wholedata from the multivar flux archive
spots=24
#this number only means the original dataset size
N = 11322
k = 32
spotchoice = rand(1:70,spots)
packages = packs(spotchoice)
varlist = ["TA_ERA"]
wholedata = reshape(hcat([Float32.(load("/net/scratch/lschulz/FLUXERAI/F32_$(flux).jld2")["data"][:,spotchoice]) for flux in varlist]...),(N,spots,length(varlist)))
wholedata=reshape(wholedata,N,spots)
outdir = "/net/scratch/lschulz/preprocessing/"
nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")
end

Threads.@threads for p=packages
    SSA_package_handler(wholedata,p,outdir,nameslist)
end

