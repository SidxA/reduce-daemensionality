"""
protophase extraction in the fluxnet file
"""

function protophase_by_varind_recm(datadir,loc_ind,season_ind,var_ind,rec_m)
    if rec_m == "raw"
        signal = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["signal"]
    elseif rec_m == "ssa"
        signal = sum(load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"],dims=2)[:,1]
    elseif rec_m == "diff"
        signal = sum(load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"],dims=2)[:,1]
    end
    return atan.(imag(hilbert_transform(Float64.(signal))),real(hilbert_transform(Float64.(signal))))
end

function locationmap(loc_ind,datadir,season_ind,rec_m,A_ind,B_ind)
    phasediff = protophase_by_varind_recm(datadir,loc_ind,season_ind,A_ind,rec_m) .- protophase_by_varind_recm(datadir,loc_ind,season_ind,B_ind,rec_m)
    h = fit(Histogram,phasediff,bins)
    h = normalize(h,mode=:pdf).weights
    return h
end


"""
phase differences between two variables at a single spot using cauchy function fit
"""

function cauchy(x,p) # fit cauchy / delta to a normalized phase difference
    x0 = p[1]
    gamma = p[2]
    @. return 1/(pi*gamma*(1+((x-x0)/gamma)^2))
end

function fit_cauchy(onebins) # perform the fit to the histogram list
    p0 = ones(2)
    return coef(curve_fit(cauchy,Array(bins)[1:end-1],onebins,p0))
end

function phase_diff(datadir,loc_ind,season_ind,rec_m,ind_A,ind_B) #loc ind needs to be a list
    bins = (-pi:0.01:pi)
    bin_N = length(bins)-1
    args = [datadir,season_ind,rec_m,ind_A,ind_B]
    binlist = map(x -> locationmap(x,args...), loc_ind)
    cauchypar = reshape(hcat(map(fit_cauchy,binlist)...),2,length(loc_ind))
    x0 = cauchypar[1,:]
    gamma = cauchypar[2,:]
    return x0,gamma
end

"""
THE INDIVIDUAL MODES and their frequencies
"""

function protofrequency(signal::Vector{Float32})
    protophase = atan.(imag(hilbert_transform(Float64.(signal))),real(hilbert_transform(Float64.(signal))))
    return Float32(2*length(signal) / count_sign_flip(protophase))
end

function mode_protofrequencies(loc_ind,datadir,season_ind,var_ind,rec_m)
    if rec_m == "ssa"
        RC = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"]
        lambda = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["lambda"]
    elseif rec_m == "diff"
        RC = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"]
        lambda = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["lambda"]
    end
    frequencies = 365 ./[protofrequency(RC[:,k]) for k in 1:48]

    return frequencies::Vector{Float32}
end

function list_harmonics(hmax)
    L = Float32[]
    for h=1:hmax
        L = append!(L,1/h)
        L = append!(L,h)
    end
    return L
end

"""
show the frequencies and phase differences for all spots
in a concise manner

"""
function flux_results(loc_ind,savename,rec_m)
    F = Figure(resolution = (1600, 600))
    datadir ="data/fluxfull/"
    season_ind = 1
    x = 1:length(loc_ind)
    axis_list = []
    """
    the individual modes
    """
    h = list_harmonics(20)
    varnames = ["T","VPD","P"]
    axis_slots =  [[1,1],[2,2],[1,3]]
    colors = [:red,:blue,:orange]
    for i = 1:3
        var_ind = i
        slot = axis_slots[i]
        varname = varnames[i]
        c=colors[i]
        f = [mode_protofrequencies(loc,datadir,season_ind,var_ind,rec_m) for loc in loc_ind]
        f = hcat([ff for ff in f]...)
        a =Axis(F[slot...],yscale = log10,ylabel=L"f [a^{-1}]",xlabel=L"locations",title = varname)
        hidespines!(a)
        a.aspect = AxisAspect(2)
        hlines!(a, h, color = :black,linewidth=1,linestyle=:dash)
        for i=1:48
            Makie.scatter!(a,x,f[i,:],markersize=8,color=c,marker=:x)
        end
        axis_list = push!(axis_list,a)
    end
    """
    the whole phase differences
    """
    varnames = ["T_VPD","T_P","VPD_P"]
    var_indeces = [(1,2),(1,3),(2,3)]
    axis_slots =  [[2,1],[1,2],[2,3]]
    bins = (-pi:0.01:pi)
    bin_N = length(bins)-1
    colors = [:purple,:darkorange,:green]
    for i in 1:3
        ind_A= var_indeces[i][1]
        ind_B = var_indeces[i][2]
        varname = varnames[i]
        slot = axis_slots[i]
        c = colors[i]
        a = Axis(F[slot...],
        xlabel = L"locations", 
        ylabel = L"\Delta \tilde{\phi}",
        title = varname)
        hidespines!(a)
        a.aspect = AxisAspect(2)
        y,ye = phase_diff(datadir,loc_ind,season_ind,rec_m,ind_A,ind_B)
        errorbars!(a,x, y, ye, color = c)
        Makie.scatter!(a,x, y, markersize = 3, color = :black)
        axis_list = push!(axis_list,a)
    end
    for i in [(1,4),(2,5),(3,6)]
        linkxaxes!(axis_list[i[1]],axis_list[i[2]])
    end
    colsize!(F.layout, 1, Auto(1.2))
    colsize!(F.layout, 2, Auto(1.2))
    colsize!(F.layout, 3, Auto(1.2))
    save(dir*savename*".png",F)
end

"""
EXECUTION
"""

using CairoMakie
using MathTeXEngine

filelist = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")[5:end-3]
nameslist = String[]
for file in filelist
    filenamesplit = split(file,"_")
    years = split(filenamesplit[end-1],"-")
    period = parse(Int,years[2])-parse(Int,years[1])
    if filenamesplit[4] == "ERAI" && filenamesplit[5] == "DD" && filenamesplit[end] =="beta-3.csv" && period > 30
    nameslist = push!(nameslist,filenamesplit[2])
    end
end

x = reshape(nameslist[loc_ind],1,length(loc_ind))

"""
sort the loc_ind by latitude or by longitude of the country!
"""
table = CSV.read("logs/KW_29/"*"fluxsites.csv",DataFrame,header=0)
loc_ind=1:71

i = [findfirst(x->x[1:2]==name[1:2],table[:,1]) for name in nameslist]
i = [isnothing(ix) ? 0 : ix for ix in i]
loc_ind = loc_ind[findall(!iszero,i)]
i = i[findall(!iszero,i)]

long_ind        = loc_ind[sortperm(table[i,6])]
long            = sort(table[i,6])
lat_ind         = loc_ind[sortperm(-table[i,5])]
lat             = sort(-table[i,5])

#do we need this?

spots_raw   = 71
N_raw       = 11322
N           = 11895
W           = 5844
k           = 48
varlist     = ["TA_ERA","VPD_ERA","PA_ERA"]
datadir     = "data/fluxfull/"
datanames   = readdir(datadir)

"""
old individual
"""

# read in by index of location and variable
# the nine methods are in 3 different files
# give protophase
# make difference for the individual var-diff into a histogram struct vector with fixed size

function histvect(varA,varB) #not complete but almost
    #fixed size for the bins
    bins = (-pi:0.01:pi)
    bin_N = length(bins)-1
    #in the end we should have a big array 3(seasons)x3(recM)x3(mean,median,variance)xbin_number
    M = Array{Float32}(undef,3,3,3,bin_N)
    #produce list of filenames
    datadir
    A_ind = 1
    B_ind = 2
    #3 seasons per 3 reconstruction methods (beeing in two files) with a struct_vect of fixed size
    #iterate the nine: map the 71 each into a struct, and give back an array
    for season_ind in 1:3,(rec_ind,rec_m) in enumerate(["raw","ssa","diff"])
        #not perfect but pretty good
        args = [datadir,season_ind,rec_m,A_ind,B_ind]
        binlist = map(x -> locationmap(x,args...), locs)
        M[season_ind,rec_ind,1,:] = Float32.([mean([binlist[i][j] for i in locs]) for j in 1:bin_N])
        M[season_ind,rec_ind,2,:] = Float32.([median([binlist[i][j] for i in locs]) for j in 1:bin_N])
        M[season_ind,rec_ind,3,:] = Float32.([var([binlist[i][j] for i in locs]) for j in 1:bin_N])
    end
    return M
end

function transform_modes(loc_ind,title)
    x = 1:length(loc_ind)
    F = Figure()
    rec_m = "diff"
    for season_ind = 1:3 , var_ind = 1:3
        f = [mode_protofrequencies(loc,datadir,season_ind,var_ind,rec_m) for loc in loc_ind]
        f = hcat([ff for ff in f]...)
        a =Axis(F[var_ind,season_ind],yscale = log10,ylabel=L"f [a^{-1}]",xlabel=L"locations",title = "s$season_ind,v$var_ind")
        h = list_harmonics(20)
        hlines!(a, h, color = :black,linewidth=1,linestyle=:dash)
        for i=1:48
            Makie.scatter!(a,x,f[i,:],markersize=14,color=:red,marker=:x)
        end
    end
    save(dir*title*".png",F)
end

function makieplot(f,loc_ind,title,xvalues,season_ind,rec_ind,rec_m,ind_A,ind_B)
    Axis(f[season_ind, rec_ind],
    xlabel = "spot number", 
    ylabel = L"\Delta \tilde{\phi}",
    title = title
    )
    x = 1:length(loc_ind)
    y,ye = phase_diff(datadir,loc_ind,season_ind,rec_m,ind_A,ind_B)
    errorbars!(x, y, ye, color = :red)
    Makie.scatter!(x, y, markersize = 3, color = :black)
end

function plotall(loc_ind,xvalues,sortingname)
    varnames = ["T_VPD","T_P","VPD_PA"]
    f = Figure()
    for (season_ind,season_name) = enumerate(["lSSA","gSSA","nodeseason"]),
        (rec_ind,rec_m) in enumerate(["raw","ssa","diff"]),
        (varnamenum,var_ind) in enumerate([(1,2),(1,3),(2,3)])
                ind_A = var_ind[1]
                ind_B = var_ind[2]
                title = "$(season_name)"#_$(rec_m)_$(varnames[varnamenum])"
                makieplot(f,loc_ind,title,xvalues,season_ind,rec_ind,rec_m,ind_A,ind_B)
    end
    Makie.save(dir*sortingname*".png",f)
end