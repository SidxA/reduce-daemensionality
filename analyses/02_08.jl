include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

spot = 19
outdir="/net/scratch/lschulz/fluxfullset/"
#meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 1
L = 5478
bins = 10 .^Array(-2:0.01:2)
# Number of points
N= 5478 
# Sample period
Ts = 1 / 365
# Start time 
t0 = 0 
tmax = t0 + (N-1) * Ts
# time coordinate
t = t0:Ts:tmax



"""
the fft spectrum
"""

F = Figure()
xtic=1:6
ax_hist = Axis(F[1,1],xscale=log10,ylabel=L"\sum \lambda",xlabel="f [a⁻¹]",yscale=log10, limits=(1/9, 10^1.5,10^-4,0.9),xticks=xtic)
ax_freq = Axis(F[1,2], limits=(1/10, 10,1,10000),yscale=log10,xticks=xtic)

signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]


Four = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

lines!(ax_freq,freqs[2740:end],abs.(Four)[2740:end],color="grey")

Four ./= maximum(abs.(Four))
lines!(ax_hist,freqs[2740:end],abs.(Four)[2740:end],color="grey")
save(dir*"test.png",F)

"""
the metadata
"""
meta = CSV.read("/net/home/lschulz/logs/KW_38/meta_sites_unique.csv",DataFrame)

datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"
nameslist = readdir(datadir)
i=0
L = String[]
for name in nameslist[5:end-2]
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    metaindex = findall(x->x==String(filenamesplit[2]),meta[!,"SITE_ID"])
    years = parse(Int,filenamesplit[end-1][end-3:end]) - parse(Int,filenamesplit[end-1][end-8:end-5])
    y = years >= 15
    if  f && d && y && !isempty(metaindex)
        L = push!(L,name)
        i+=1
    end
end

indices = [4,5,13,18,24]
L_i = L[deleteat!(Array(1:25),indices)]
variable_list = ["GPP_DT_VUT_MEAN","NEE_VUT_MEAN","RECO_DT_VUT_MEAN","TS_F_MDS_1","TA_F","VPD_F"] #"SWC_F_MDS_1"
datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"
savedir = "/net/scratch/lschulz/fluxnetfullset/"

datafile = Array{Float32}(undef,5478,length(L_i),length(variable_list))
metaindices = Int64[]
for (l,name) in enumerate(L_i)
    csv = CSV.read(datadir*name,DataFrame)
    for (i,v) in enumerate(variable_list)
        datafile[:,l,i] = csv[end-5478+1:end,v]
    end
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    metaindices = push!(metaindices,findall(x->x==String(filenamesplit[2]),meta[!,"SITE_ID"])[1])
    #println(meta[metaindex,:])
    
end

sitesinfo = DataFrame(meta[metaindices,:])

jldsave(dir*"meta.jld2",meta = sitesinfo)
