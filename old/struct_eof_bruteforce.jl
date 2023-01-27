using Distributed

@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")
    

    mutable struct EOF_Diff#N,W,k
        #the variables
        emb::Matrix{Float64}
        W::Int64
        k::Int64
        DM::DiffMap{Float64}
        EOF::Matrix{Float64}
        #the init function
        function EOF_Diff(N::Int,W::Int,k::Int)
            emb             = Array{Float64}(undef,N-W+1,W) #P x W
            EOF             = Array{Float64}(undef,W,k)
            return new(emb,W,k,fit(DiffMap,zeros(2,2)),EOF)
        end
    end

    #with boundary creation effect
    #fixed W and k, list of loc (usually one)
    function Diff_package_handler(data::Matrix{Float64},package,outdir::String,names::Vector{String})
        W = package[1] # needs to be fixed
        k = package[2] # needs to be fixed
        loc_ind = package[3] # needs to be a list
        aet_combinations = package[4]
        object = EOF_Diff(12052,W,k)
        for i in loc_ind
            name = names[i]
            signal = boundary(data[:,i],365) |> centralizer
            object.emb = embed_lag(signal,object.W)'

            for aet in aet_combinations
                object.EOF = Matrix(ManifoldLearning.transform(fit(DiffMap,object.emb::Matrix{Float64},maxoutdim=object.k,t=aet[3], α=aet[1], ɛ=aet[2]))')
                jldsave(outdir*"eof_Diff_$(i)_$(name)_$(W)_$(k)_$(aet[1])_$(aet[2])_$(aet[3])_$(rand(1000:9999)).jld2",data = Matrix(object.EOF))
            end
        end
        
    end

    #create a number of packages for each name 
    function Diff_package_creation(W::Int64,k::Int64,aet::Vector{Tuple{Int64,Int64,Int64}},names_ind::Vector{Int64},N_p::Int64)
        L = length(names_ind)
        packages = []

        for p in 1:L
            for j = 1:N_p
                packages = push!(packages,(W,k,p,aet[j:N_p:end]))
            end
        end
        return packages
    end


    wholedata = Array(CSV.read("/net/scratch/lschulz/data_deseason",types=Float64,DataFrame))
    nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")
    outdir =  "/net/scratch/lschulz/diff_ttest/"

    k = 32
    W = 4383
    aet =[(1,64,1),(1,64,4),(1,64,10),(1,64,20),(1,64,50),(1,64,100)]
    names_ind = Array(1:length(nameslist))
    N_p = 2

    packages = Diff_package_creation(W,k,aet,names_ind,N_p)

end

Threads.@threads for p=packages
    Diff_package_handler(wholedata,p,outdir,nameslist)
end

