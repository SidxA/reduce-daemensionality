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

tw = t0:Ts:(t0+(W-1)*Ts)



"""
the individual modes and their fft spectrum
"""
meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]
savedir=dir*"spectra/"


outdir="/net/scratch/lschulz/fluxfullset_1a/"
W = 365
freqstart = Int(round(W/2,digits=0))+1
tw = t0:Ts:(t0+(W-1)*Ts)

for spot = 1:20
name = "$(spot)_1a"
F = Figure(resolution=(800,800))

xtic=1:6
title1="$spot  $(meta[spot,"name"]) $(meta[spot,"IGBP_class"]) W = $(round(W/365,digits=0)) a"
subtitle1 = "time series"

#axis

ax_signal = Axis(F[1,1:2],xlabel="time",ylabel="GPP",yticksvisible = false,
yticklabelsvisible = false,title = title1,
subtitle = subtitle1,
titlealign = :left,
subtitlegap = 2,
titlegap = 5,
subtitlefont = :italic,
subtitlelineheight = 0.9,
titlelineheight = 0.9,)

ax_modes = Axis(F[2,1],xlabel="time",ylabel="GPP",yticksvisible = false,
yticklabelsvisible = false,title = "SSA modes",
titlealign = :left,
subtitlegap = 2,
titlegap = 5,
subtitlefont = :italic,
subtitlelineheight = 0.9,
titlelineheight = 0.9,)


ax_freq = Axis(F[2,2],limits=(1/10, 7,0.001,17),#,yscale=log10
xticks=xtic,xlabel="frequency",ylabel ="rel power",yticksvisible = false,
yticklabelsvisible = false,title = "SSA spectrum",
titlealign = :left,
subtitlegap = 2,
titlegap = 5,
subtitlefont = :italic,
subtitlelineheight = 0.9,
titlelineheight = 0.9,)

ax_modes2 = Axis(F[3,1],xlabel="time",ylabel="GPP",yticksvisible = false,
yticklabelsvisible = false,title = "NLSA modes",
titlealign = :left,
subtitlegap = 2,
titlegap = 5,
subtitlefont = :italic,
subtitlelineheight = 0.9,
titlelineheight = 0.9,)


ax_freq2 = Axis(F[3,2],limits=(1/10, 7,0.001,17),#,yscale=log10
xticks=xtic,xlabel="frequency",ylabel ="rel power",yticksvisible = false,
yticklabelsvisible = false,title = "NLSA spectrum",
titlealign = :left,
subtitlegap = 2,
titlegap = 5,
subtitlefont = :italic,
subtitlelineheight = 0.9,
titlelineheight = 0.9,)

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

    lines!(ax_modes,1:W,mode .+ (k*0.07))
    lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k)
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

    lines!(ax_modes2,1:W,mode .+ (k*0.07))
    lines!(ax_freq2,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k)
end

colgap!(F.layout,1,0)
rowgap!(F.layout,1,0)
rowgap!(F.layout,2,0)

save(savedir*name*".png",F)

end
