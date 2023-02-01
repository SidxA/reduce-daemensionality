
include("/net/home/lschulz/scripts/toolbox.jl")
include("/net/home/lschulz/scripts/struct_blockdata.jl")
include("/net/home/lschulz/scripts/lSSA_deseasonalization.jl")

#big parameters
#N = 11322
#W = Int(floor((365.25*11)))
k = 48
"""
running raw at the moment
"""
#deseasonalization
#windowlength = 30
#overlap=10

#fluxnetfullset
wholedata = SharedArray{Float32}(load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts_ta_vpd.jld2")["wholedata"])
#LAIF
#wholedata = SharedArray{Float32}(Matrix(load("/net/scratch/lschulz/LAIF/data.jld2")["data"]'))
N,spots,nvars = size(wholedata)
vars=1

"""
OUTDIR NEEDS TO HAVE /
"""

outdir = "/net/scratch/lschulz/fluxfullset_1a/"
#outdir="/net/scratch/lschulz/LAIF/"

#spots=71
#varlist = ["TA_ERA"]#,"VPD_ERA","PA_ERA"]
#wholedata = reshape(hcat([Float32.(load("/net/scratch/lschulz/FLUXERAI/F32_$(flux).jld2")["data"]) for flux in varlist]...),(N,spots,length(varlist)))
#outdir = "/net/scratch/lschulz/ta_erai_11a/"

"""
put EVERYTHING variable TO KNOW in the package
"""
#loc_ind var_ind
struct package
    outdir::String
    loc_ind::Int64
    var_ind::Int64
    function package(outdir,loc_ind,var_ind)
        return new(outdir,loc_ind,var_ind)
    end
end

function packaging(outdir)
    if isdir(outdir)==false
        mkdir(outdir)
    end
    P = package[]
    for spot in 1:spots,vari in 1:vars
        P = push!(P,package(outdir,spot,vari))
    end
    return P
end



#iterate the packages of all specified parameters

function doitall(d::matrixholder,wholedata::SharedArray{Float32},p::package)
    outdir = p.outdir
    loc_ind= p.loc_ind
    var_ind= p.var_ind

    #preproc year manual
    #sampleyear = 24
    sampleyear = 365

    put_data_raw(d,wholedata[:,loc_ind,var_ind])

    put_epsilon(d)
    put_EOF_diff(d)
    calculate(d)
    save_results(d,outdir*"diff_$(d.W)_$(loc_ind)_$(var_ind)_raw")

    put_EOF_ssa(d)
    calculate(d)
    save_results(d,outdir*"ssa_$(d.W)_$(loc_ind)_$(var_ind)_raw")

    #put_data_gSSA(d,wholedata[:,loc_ind,var_ind]::Vector{Float32},sampleyear)

    #put_epsilon(d)
    #put_EOF_diff(d)
    #calculate(d)
    #save_results(d,outdir*"diff_$(d.W)_$(loc_ind)_$(var_ind)_gSSA")

    #put_EOF_ssa(d)
    #calculate(d)
    #save_results(d,outdir*"ssa_$(d.W)_$(loc_ind)_$(var_ind)_gSSA")

    #put_data_lSSA(d,wholedata[:,loc_ind,var_ind]::Vector{Float32},sampleyear)

    #put_epsilon(d)
    #put_EOF_diff(d)
    #calculate(d)
    ##save_results(d,outdir*"diff_$(d.W)_$(loc_ind)_$(var_ind)_lSSA")

    #put_EOF_ssa(d)
    #calculate(d)
    #save_results(d,outdir*"ssa_$(d.W)_$(loc_ind)_$(var_ind)_lSSA")

    #put_data_window(d,wholedata[:,loc_ind,var_ind]::Vector{Float32},sampleyear)

    #put_epsilon(d)
    #put_EOF_diff(d)
    #calculate(d)
    #save_results(d,outdir*"diff_$(d.W)_$(loc_ind)_$(var_ind)_win")

    #put_EOF_ssa(d)
    #calculate(d)
    #save_results(d,outdir*"ssa_$(d.W)_$(loc_ind)_$(var_ind)_win")

end

#need to manually put the sampleyear in here...
#Walist = [5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7]
#Walist = [15,16,17]

Threads.@threads for p=packaging(outdir)
    for W = 365#Int.(floor.(Walist .* 24))
        d = matrixholder(N,W,k)
        doitall(d,wholedata,p)
    end
end
