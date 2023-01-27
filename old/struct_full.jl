#this is the package based approach with 1 struct for both algorithms
#the struct holds the embed lag, eof, pc and rc
# and then it saves it to one file

# N, W and k need to be specified

#without boundary creation effect - simple reconstruction!
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

    #fixed N, W and k, list of loc
    #then iterates the parameters
    #with centralization!
    #saves as name number N W k 

    function Diff_package_handler(data::Matrix{Float32},package,outdir::String,nameslist::Vector{String})
        N = package[1]
        W = package[2] # needs to be fixed
        k = package[3] # needs to be fixed
        P = N - W +1
        loc_ind = package[4] # needs to be a list
        aet_combinations = package[5]
        object = matrix_holder(N,W,k,data[:,loc_ind])
        for i in loc_ind
            name = nameslist[i]

            object.emb = centralized_embed_lag(object.signal,W)'

            for aet in aet_combinations
                object.EOF = Matrix(ManifoldLearning.transform(fit(DiffMap,object.emb::Matrix{Float32},maxoutdim=object.k,t=aet[3], α=aet[1], ɛ=aet[2]))')
                object.PC = hcat([pc!(object.signal,object.EOF[:,i],P,W) for i in 1:k]...)'
                object.RC = hcat([reconstructor(object.PC[i,:],object.EOF[:,i],N,W) for i in 1:k]...)
                jldsave(outdir*"Diff_$(i)_$(name)_$(N)_$(W)_$(k)_$(aet[1])_$(aet[2])_$(aet[3])_$(rand(1000:9999)).jld2",
                EOF = Matrix(object.EOF),PC = Matrix(object.PC),RC = Matrix(object.RC))
            end
        end
        
    end

    function SSA_package_handler(data::Matrix{Float32},package,outdir::String,nameslist::Vector{String})
        N = package[1]
        W = package[2] # needs to be fixed
        k = package[3] # needs to be fixed
        P = N - W +1
        loc_ind = package[4] # needs to be a list
        aet_combinations = package[5]
        object = matrix_holder(N,W,k,data[:,loc_ind])
        for i in loc_ind
            name = nameslist[i]

            object.emb = centralized_embed_lag(object.signal,W)'
            object.EOF = Matrix(projection(fit(PCA,object.emb',method=:svd,pratio=1.0,maxoutdim=k)))
            object.PC = hcat([pc!(object.signal,object.EOF[:,i],P,W) for i in 1:k]...)'
            object.RC = hcat([reconstructor(object.PC[i,:],object.EOF[:,i],N,W) for i in 1:k]...)

            jldsave(outdir*"SSA_$(i)_$(name)_$(N)_$(W)_$(k)_$(0.0)_$(0.0)_$(0.0)_$(rand(1000:9999)).jld2",
            EOF = Matrix(object.EOF),PC = Matrix(object.PC),RC = Matrix(object.RC))

        end
        
    end
    #create a number of packages for each name
    function Diff_package_creation(N::Int64,W::Int64,k::Int64,aet::Vector{Tuple{Int64,Int64,Int64}},names_ind::Vector{Int64},N_p::Int64)
        L = length(names_ind)
        packages = []

        for p in 1:L
            for j = 1:N_p
                packages = push!(packages,(N,W,k,p,aet[j:N_p:end]))
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

    wholedata = Float32.(Array(CSV.read("/net/scratch/lschulz/data_deseason",types=Float64,DataFrame)))
    nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")
    outdir =  "/net/scratch/lschulz/ssasweep/" #always put the / at the end!

    N = 12052
    k = 32
    W = 1461
    #aet =[(1,64,1),(1,64,4),(1,64,10),(1,64,20),(1,64,50),(1,64,100)]
    aet=[(0,0,0)]
    names_ind = Array(1:length(nameslist))
    N_p = 1

    function create_p_by_W()
        #packages = Diff_package_creation(N,W,k,aet,names_ind,N_p)
        running_packages = []

        for W in Int.(floor.(Array(365.25:365.25:365.25*10)))
            running_packages= append!(running_packages,Diff_package_creation(N,W,k,aet,names_ind,N_p))
        end
        return running_packages
    end
    #always need to give the parameter packages as a function instance, not variable!
    running_packages=create_p_by_W()
end

Threads.@threads for p=running_packages
    SSA_package_handler(wholedata,p,outdir,nameslist)
end

