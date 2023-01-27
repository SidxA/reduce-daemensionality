
using Distributed


@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")
    

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

    function SSA_package_handler(data::Array{Float32},package)
        N = package[1]
        W = package[2] # needs to be fixed
        k = package[3] # needs to be fixed
        P = N - W +1
        loc_index = package[4]
        aet_combinations = package[5]
        var_index = package[6]
        outdirname = package[7]
        spotchoice_index = package[8]

        object = matrix_holder(N,W,k,data[:,spotchoice_index,var_index])


        object.emb = centralized_embed_lag(object.signal,W)'
        object.EOF = Matrix(projection(fit(PCA,object.emb',method=:svd,pratio=1.0,maxoutdim=k)))
        object.PC = hcat([pc!(object.signal,object.EOF[:,i],P,W) for i in 1:k]...)'
        object.RC = hcat([reconstructor(object.PC[i,:],object.EOF[:,i],N,W) for i in 1:k]...)

        jldsave(outdirname,
        EOF = Matrix(object.EOF),PC = Matrix(object.PC),RC = Matrix(object.RC))

        
    end

    #create a number of packages for each name
    function package_creation(N::Int64,W::Int64,k::Int64,spotchoice::Vector{Int64},varlist,outdir,nameslist)
        packages = []
        for (spotchoice_ind,loc_ind) in enumerate(spotchoice)
            name = nameslist[loc_ind]
            for varindex in 1:length(varlist)
                outdirname = outdir*"SSA_$(loc_ind)_$(name)_$(N)_$(W)_$(k)_$(varindex)_$(varlist[varindex]).jld2"
                packages = push!(packages,(N,W,k,loc_ind,(0,0,0),varindex,outdirname,spotchoice_ind))
            end
        end
        return packages
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

    spots=24
    N = 11322
    k = 32

    spotchoice = rand(1:70,spots)
    varlist = ["TA_ERA","SW_IN_ERA","LW_IN_ERA","VPD_ERA","PA_ERA","P_ERA","WS_ERA"]
    wholedata = reshape(hcat([Float32.(load("/net/scratch/lschulz/FLUXERAI/F32_$(flux).jld2")["data"][:,spotchoice]) for flux in varlist]...),(N,spots,length(varlist)))
    outdir = "/net/scratch/lschulz/multivar/"
    nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")

    W1 = 4383
    W2 = 5113
    W3 = 2922
    running_packages=append!(package_creation(N,W1,k,spotchoice,varlist,outdir,nameslist),package_creation(N,W2,k,spotchoice,varlist,outdir,nameslist),package_creation(N,W3,k,spotchoice,varlist,outdir,nameslist))
end



Threads.@threads for p=running_packages
    SSA_package_handler(wholedata,p)
end
