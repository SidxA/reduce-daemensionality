using Distributed

@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")
    
    #for one fixed W and N we iterate over the locations!
    mutable struct EOF_SSA
        #the variables
        emb::Matrix{Float64}
        W::Int64
        k::Int64
        pca::PCA{Float64}
        EOF::Matrix{Float64}
        #the init function
        function EOF_SSA(N::Int,W::Int,k::Int)
            emb             = Array{Float64}(undef,W,N-W+1)
            EOF             = Array{Float64}(undef,W,k)
            return new(
                emb,W,k,fit(PCA,zeros(2,2)),EOF
            )
        end
    end

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


    function perform_SSA(data::Matrix{Float64},outdir::String,names::Vector{String},W::Int,k::Int)
        object = EOF_SSA(12052,W,k)
        for (i,name) in enumerate(names)
            signal = data[:,i]
            object.emb = embed_lag(signal,object.W)::Matrix{Float64}
            object.pca = fit(PCA,object.emb,method=:svd,pratio=1.0,maxoutdim=object.k)
            object.EOF = projection(object.pca)
            jldsave(outdir*"eof_SSA_$(name)_$(W)_$(k)_0_0_0_$(rand(1000:9999)).jld2",data = Matrix(object.EOF))
        end
    end

    function perform_Diff(data::Matrix{Float64},outdir::String,names::Vector{String},W::Int,k::Int,aet_combinations)
        object = EOF_Diff(12052,W,k)
            for (i,name) in enumerate(names)
                signal = boundary(data[:,i],365) |> centralizer
                #setfield!(object,:emb,embed_lag(signal,object.W))
                object.emb = embed_lag(signal,object.W)'
                for aet in aet_combinations
                    #object.DM = fit(DiffMap,object.emb::Matrix{Float64},maxoutdim=object.k,t=aet[3], α=aet[1], ɛ=aet[2])
                    #object.EOF = Matrix(ManifoldLearning.transform(object.DM)'){Float64}
                    object.EOF = Matrix(ManifoldLearning.transform(fit(DiffMap,object.emb::Matrix{Float64},maxoutdim=object.k,t=aet[3], α=aet[1], ɛ=aet[2]))')

                    jldsave(outdir*"eof_Diff_$(name)_$(W)_$(k)_$(aet[1])_$(aet[2])_$(aet[3])_$(rand(1000:9999)).jld2",data = Matrix(object.EOF))
                end
            end
    end


    function parameter_wrapper(wholedata,nameslist,N_W,Wlist,N_loc) #[(SSA/Diff,W,names,0/aet)]
        aetlist = [(1.0,64,1),(1.0,64,2),(1.0,64,4),(1.0,64,7),(1.0,64,10),(1.0,64,15),(1.0,64,20)]
        #aetlist = [(1.0,64,1),(1.0,128,1),(1/2,64,1),(1/2,128,1),(1.0,64,5),(1.0,128,5),(1/2,64,5),(1/2,128,5)]#,(1.0,64,10),(1.0,128,10),(1/2,64,10),(1/2,128,10)]
        para = []
        L = length(nameslist)
        for word in ["Diff"]#,"SSA"]
            for W in rand(Wlist,N_W)
                indices = rand(1:L,N_loc)
                randnames = nameslist[indices]
                data = wholedata[:,indices]
                para = push!(para,(word,data,randnames,W,aetlist))
            end
        end
        return para
    end

    wholedata = Array(CSV.read("/net/scratch/lschulz/data_deseason",types=Float64,DataFrame))
    nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")
    outdir =  "/net/scratch/lschulz/diff_t/"
    #Wlist = Int64[Int.(floor.((365.25*5):91.3125:6026).+1)...]
    #Wlist = Int64[Int.(floor.((365.25*12):91.3125:6026).+1)...]
    Wlist = [730,2922,3652]

    k = 32
    parameterlist = parameter_wrapper(wholedata,nameslist,3,Wlist,20)
end

Threads.@threads for p=parameterlist
    method = p[1]
    data = p[2]
    names = p[3]
    W = p[4]
    aetlist = p[5]
    if method == "SSA"
        perform_SSA(data,outdir,names,W,k)
    elseif method == "Diff"
        perform_Diff(data,outdir,names,W,k,aetlist)
    elseif time()>= 1.660289069097841e9
        break
    end
end


    #=

    #for one W and k and nameslist open file, init struct, loop over locations, close
    function perform_SSA(outdir::String,data,names,W::Int,k::Int)
        #outfile=outdir*"eofs_ssa_$(W)_$(k)_$(rand(1000:9999)).jld2"
        object = EOF_SSA(12052,W,k)
        #f = jldopen(outfile, "a+")::JLD2.JLDFile{JLD2.MmapIO}
        #jldopen(outfile, "a+") do f
            for name in names
                # do something
                object.emb = embed_lag(data[!,name],object.W)::Matrix{Float64}
                object.pca = fit(PCA,object.emb,method=:svd,pratio=1.0,maxoutdim=object.k)
                object.EOF = projection(object.pca)
                jldsave(outdir*"eof_S_$(name)_$(W)_$(k)_0_0_0_$(rand(1000:9999)).jld2",Matrix(object.EOF))

                #write(f,name*"_$(rand(1000:9999))",object.EOF::Matrix{Float64})
            end
        #end
        #close(f)

    end

    function perform_Diff(outdir::String,data::DataFrame,names,W::Int,k::Int,aet_combinations)
        #outfile = outdir*"eof_D_$(W)_$(k)_$(rand(1000:9999)).jld2"
        object = EOF_Diff(12052,W,k)

        #jldopen(outfile, "a+") do f
        #f = jldopen(outfile, "a+")::JLD2.JLDFile{JLD2.MmapIO}
            for name in names
                # do something
                object.emb = embed_lag(data.name::Vector{Float64},object.W)::Matrix{Float64}
                for aet in aet_combinations

                    object.DM = fit(DiffMap,object.emb,maxoutdim=object.k,t=aet[3], α=aet[1], ɛ=aet[2])
                    object.EOF = Matrix(ManifoldLearning.transform(object.DM)')
                    #write(f,name*string("_",aet)*"_$(rand(1000:9999))",object.EOF::Matrix{Float64})
                    jldsave(outdir*"eof_D_$(name)_$(W)_$(k)_$(aet[1])_$(aet[2])_$(aet[3])_$(rand(1000:9999)).jld2",Matrix(object.EOF))
                end
            end
            #close(f)
        #end
    end
=#