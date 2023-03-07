
include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

"""
recreate the bigf figure
"""
#uses this function
function single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples)
    N = 5478
    k = 48
    filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    lambda,protofreq,RC = extract_from_single_file(filename,yearsamples,N,k)
    lambda_l = sort(lambda,rev=true)[1:kappa]
    protofreq_l = protofreq[sortperm(lambda,rev=true)][1:kappa]
    RC_l = RC[:,sortperm(lambda,rev=true)][:,1:kappa]
    return RC_l,lambda_l,protofreq_l
end

spot = 3
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
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
Ts = 1 / 365.25
# Start time 
t0 = 0 
tmax = t0 + (N-1) * Ts
# time coordinate
t = t0:Ts:tmax

xtic = 1:7#round.(Array([1/7:1/7:5/7...,[1,2,3,4,7] ...]),digits=2)

color_ssa = "darkgreen"
color_nlsa = "purple"
ms = 18

F = Figure(resolution= (800,500))
ax_fft = Axis(F[1,1:2],limits=(1/9,7,0.0001,1),xticks=xtic,yscale=log10)
ax_both = Axis(F[2,1:2],yscale=log10,#xscale=log10,
limits=(1/9,7,10^-4,0.9),xticks=xtic,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(7),)
ax_spec = Axis(F[2,3]
,yscale=log10,limits=(0,49,10^-4,0.9),)

#colsize!(F.layout, 1, Fixed(250))
#colsize!(F.layout, 2, Fixed(250))
colsize!(F.layout, 3, Fixed(120))

#rowsize!(F.layout, 1, Fixed(100))
#rowsize!(F.layout, 2, Fixed(200))

harm = list_pure_harmonics(7)

vlines!(ax_both,harm[end-6:end],color="grey")
#vlines!(ax_hist,harm,color="grey")


signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]
Four = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
Four ./= maximum(abs.(Four))
lines!(ax_fft,freqs[2740:end], abs.(Four)[2740:end], color=:black,label="signal")


for (i,method) = enumerate(["ssa","diff"])


    r_GPP = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples) #RC,lambda,freq
    RC = r_GPP[1]
    lambda = r_GPP[2]
    protofreq = r_GPP[3]
    scatter!(ax_both,protofreq,lambda,marker=[:cross,:xcross][i],
    markersize=ms,markerstrokewidth=4,color=[color_ssa,color_nlsa][i])

    h_weights =fit(Histogram,protofreq,weights(lambda),bins).weights
    #lines!(ax_hist,bins[1:end-1],10^-5 .+ h_weights,label=["SSA","NLSA"][i])

    lines!(ax_spec,1:48,lambda,color=[color_ssa,color_nlsa][i])

    Four = fft(sum(RC,dims=2)[:]) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    Four ./= maximum(abs.(Four))
    lines!(ax_fft,freqs[2740:end], abs.(Four)[2740:end],
    color=[color_ssa,color_nlsa][i],label=["SSA","NLSA"][i])
end


hidespines!(ax_both, :t, :r)

Legend(F[1,3],ax_fft,"Method")


linkyaxes!(ax_both,ax_spec)
hidespines!(ax_fft, :t, :r, :b)
hidespines!(ax_spec, :t, :r)
hideydecorations!(ax_spec,grid=false)
linkxaxes!(ax_fft,ax_both)
linkyaxes!(ax_fft,ax_both)

hidexdecorations!(ax_fft)

colgap!(F.layout,1,0)
rowgap!(F.layout,1,6)


save(dir*"bigf.png",F)



"""
the individual modes and their fft spectrum
"""

meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]
outdir="/net/scratch/lschulz/fluxfullset/"
W = 2556
freqstart = Int(round(W/2,digits=0))+1
tw = t0:Ts:(t0+(W-1)*Ts)
spot=3

years = (1:W) ./ 365

F = Figure(resolution=(800,600))

color_ssa = "darkgreen"
color_nlsa = "purple"

xtic=1:6

#axis

ax_modes = Axis(F[1,1:2],yticksvisible = false,
yticklabelsvisible = false,xticks=1:7)


ax_freq = Axis(F[1,3:4],limits=(1/10, 7,-0.7,37),#,yscale=log10
xticks=xtic,yticksvisible = false,
yticklabelsvisible = false)


#signal

signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]


Four = fft(signal) |> fftshift
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
Four ./= maximum(abs.(Four))

lines!(ax_signal,1:N,signal,color="grey")
#lines!(ax_freq,freqs[2740:end],abs.(Four)[2740:end],color="grey")

#modes

method = "ssa"

Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]

for k = 1:16

    mode = Eof[:,k] 
    Four = fft(mode) |> fftshift
    freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
    Four ./= maximum(abs.(Four)) 

    lines!(ax_modes,years,mode .+ (k*0.07),
    color=color_ssa)
    lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k,
    color=color_ssa)
end

method = "diff"

Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]

for k = 1:16

    mode = Eof[:,k] 
    Four = fft(mode) |> fftshift
    freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
    Four ./= maximum(abs.(Four)) 

    lines!(ax_modes,years,mode .+ ((18+k)*0.07),
    color=color_nlsa)
    lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k .+ 18,
    color=color_nlsa)
end

colgap!(F.layout,1,0)

hidespines!(ax_freq, :t, :r, :l)
hidespines!(ax_modes, :t, :r, :l)

save(dir*"modeshape.png",F)

"""
plot individual signal + reconstructions
"""

method = "ssa"
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
rc_ssa = sum(file["RC"],dims=2)[:]
method = "diff"
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
rc_nlsa = sum(file["RC"],dims=2)[:]
signal = file["signal"]

color_ssa = "darkgreen"
color_nlsa = "purple"
lw = 2
ms = 5

years = ((1:5478) ./ 365) .+ 1985

F = Figure(resolution=(800,300))

ax = Axis(F[1,1],xticks=1985:3:2000,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(3),)

scatter!(ax,years,signal,color="black",markersize = ms,marker=:x)
lines!(ax,years,signal,color="black",markersize = 1)
lines!(ax,years,rc_ssa,color=color_ssa,linewidth=lw)
lines!(ax,years,rc_nlsa,color=color_nlsa,linewidth=lw)

hideydecorations!(ax)
hidespines!(ax,:t,:r)

save(dir*"reconstructions.png",F)
