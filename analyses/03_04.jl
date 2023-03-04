
include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie


"""
plot the seasonality behavior...
"""
function gauss(x,p) # cauchy peak with two parameters
    x0 = p[1]
    gamma = p[2]
    sigma = p[3]
    #1/(sigma * sqrt(2*pi))
    @. return gamma * exp.(-(x-x0)^2/sigma^2/2)
end

function fit_gauss(yvalues,xvalues) # fit the cauchy peak to a vec of numbers
    onebins = xvalues
    bins = yvalues
    p0 = ones(3) .* 5
    return coef(curve_fit(gauss,bins,onebins,p0))
end


function mode_harmonicity_estimation_gauss(Eof,threshold)
    li_harmonics = Int64[]
    li_mixed = Int64[]
    li_h_freq = Float64[]
    li_m_freq = Float64[]
    li_residual=Int64[]
    W,k = size(Eof)

    for i in 1:k

        mode = Eof[:,i]

        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift

        freqstart = findall(x->x>=1/12,freqs)[1]
        freqend = findall(x->x>=6,freqs)[1]

        
        val_domain = abs.(Four[freqstart:freqend])
        #val_domain ./= maximum(val_domain)

        freq_domain = freqs[freqstart:freqend]

        freq, value,sigma = fit_gauss(freq_domain,val_domain)
        peak = gauss(freq_domain,(freq,value,sigma))
        residual = val_domain .- peak

        freq_round = round(freq,digits=1)
        
        val = abs(1/sigma /sqrt(2*pi))
        
        #println("$i \t $freq \t $val \t $(maximum(residual))")

        #if maximal residual value other then single peak is larger then relative peak threshold value, its 
        if maximum(residual .+ 0.0)*threshold <= value &&  freq_round in Float64.(1:10)
            li_harmonics = append!(li_harmonics,i)
            li_h_freq = append!(li_h_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= value && freq_round in Float64.(1:10)
            li_mixed = append!(li_mixed,i)
            li_m_freq = append!(li_m_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= value
            li_residual = append!(li_residual,i)
        else
            #println("no peak")
        end

    end

    return li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual
end


# vielsalm: 3


threshold = 10
savedir=dir*"rec/"
variable_list = ["GPP","NEE","RECO","soil water content","soil temperature"]
meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]
outdir="/net/scratch/lschulz/fluxfullset/"

W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 1
# Number of points
N= 5478 
# Sample period
Ts = 1 / 365
# Start time 
t0 = 0 
tmax = t0 + (N-1) * Ts
# time coordinate
t = t0:Ts:tmax
0:Ts:(t0+(W-1)*Ts)
tw = t0:Ts:(t0+(W-1)*Ts)


savedirname = dir*"poster/"
spot = 3
threshold = 10
vari=1

#ssa

method="ssa"
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]

# ssa reconstruction: freq 1

li_harmonics=[1,2]
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
ssa_trend_first = sum(RC[:,li_harmonics],dims=2)[:]

# ssa reconstruction: freq harmonics

li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
ssa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

#nlsa 

method="diff"
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]

# nlsa reconstruction: freq 1

li_harmonics=[1,2]
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
nlsa_trend_first = sum(RC[:,li_harmonics],dims=2)[:]

# nlsa reconstruction: freq harmonics

li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
nlsa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

trajectorymatrix = hcat(
ssa_trend_first .+ 0.0,
ssa_trend_harm .+ 2.0,
nlsa_trend_first.+ 4.0,
nlsa_trend_harm .+ 6.0)

"""
plotting procedure : fundamental versus harmonic
"""

#signals

years = 1985 .+ ((1:5478) ./ 365)
fundamental = ssa_trend_first
harmonic = ssa_trend_harm
res_fund = signal .- fundamental
res_harm = signal .- harmonic

F = Figure(resolution=(800,600))
ax = Axis(F[1,1],
xticks = 1985:3:2000,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(3),
)

lw = 3
ms = 4
col_fundamental = "red"
col_harmonic = "blue"
col_signal = "black"

#signal background

scatter!(ax,years,signal .+ 8.0,linewidth=lw,
color=col_signal,marker=:x,markersize=ms)
scatter!(ax,years,signal .+ 12.0,linewidth=lw,
color=col_signal,marker=:x,markersize=ms,label="signal")

#fundamental
lines!(ax,years,fundamental .+ 12.0,linewidth=lw,
color=col_fundamental,linestyle=:solid,label="fundamental")
lines!(ax,years,res_fund .+ 4.0,linewidth=lw,
color=col_fundamental,linestyle=:solid,label="fundamental")
#harmonic
lines!(ax,years,harmonic .+ 8.0,linewidth=lw,
color=col_harmonic,linestyle=:solid,label="harmonic")
lines!(ax,years,res_harm .+ 0.0,linewidth=lw,
color=col_harmonic,linestyle=:solid,label="harmonic")

hideydecorations!(ax)
hidespines!(ax, :t, :r, :l) # only top and right
#axislegend(ax)

save(dir*"seasons.png",F)

#now the additional FFT

ms = 8
normalizer(x) = x./maximum(x)
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

spec_fun = (abs.(fft(fundamental) |> fftshift)[freqstart:freqend] |> normalizer ) .* 10^12
spec_res_fun = (abs.(fft(res_fund) |> fftshift)[freqstart:freqend] |> normalizer ) .* 10^4
spec_har = (abs.(fft(harmonic) |> fftshift)[freqstart:freqend] |> normalizer ) .* 10^8
spec_res_har = (abs.(fft(res_harm) |> fftshift)[freqstart:freqend] |> normalizer )
spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )

F = Figure(resolution=(800,600))
ax = Axis(F[1,1],
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(7),
yscale = log10)

#signal
scatter!(ax,freq_domain,spec_signal .*10^8,linewidth=lw,
color=col_signal,marker=:x,markersize=ms)
scatter!(ax,freq_domain,spec_signal .*10^12,linewidth=lw,
color=col_signal,marker=:x,markersize=ms,label="signal")

#fundamental
lines!(ax,freq_domain,spec_fun,linewidth=lw,
color=col_fundamental,linestyle=:solid)
lines!(ax,freq_domain,spec_res_fun,linewidth=lw,
color=col_fundamental,linestyle=:solid)

#harmonic
lines!(ax,freq_domain,spec_har,linewidth=lw,
color=col_harmonic,linestyle=:solid)
lines!(ax,freq_domain,spec_res_har,linewidth=lw,
color=col_harmonic,linestyle=:solid)

hideydecorations!(ax)
hidespines!(ax, :t, :r, :l) # only top and right
#axislegend(ax)

save(dir*"seasons_spec.png",F)


"""
plotting procedure : nlsa versus ssa
"""


F = Figure(resolution=(800,600))
ax = Axis(F[1,1],
xticks = 1985:3:2000,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(3),
)

lw = 3
col_ssa = "darkgreen"
col_nlsa = "purple"
col_signal = "black"

ms = 4

#signal background

scatter!(ax,years,signal .+ 8.0,linewidth=lw,
color=col_signal,marker=:x,markersize=ms)
scatter!(ax,years,signal .+ 12.0,linewidth=lw,
color=col_signal,marker=:x,markersize=ms,label="signal")

#ssa
lines!(ax,years,ssa_trend_harm .+ 12,linewidth=lw,
color=col_ssa,linestyle=:solid,label="ssa")
lines!(ax,years,signal .- ssa_trend_harm .+ 4,linewidth=lw,
color=col_ssa,linestyle=:solid,label="ssa")

#nlsa
lines!(ax,years,nlsa_trend_harm .+ 8,linewidth=lw,
color=col_nlsa,linestyle=:solid,label="nlsa")
lines!(ax,years,signal .- nlsa_trend_harm,linewidth=lw,
color=col_nlsa,linestyle=:solid,label="nlsa")

hideydecorations!(ax)
hidespines!(ax, :t, :r, :l) # only top and right
#axislegend(ax)

save(dir*"comp.png",F)

#now the additional FFT

ms = 8
normalizer(x) = x./maximum(x)
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

spec_ssa = (abs.(fft(ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^12
spec_res_ssa = (abs.(fft(signal .- ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^4
spec_nlsa = (abs.(fft(nlsa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^8
spec_res_nlsa = (abs.(fft(signal .- nlsa_trend_harm) |> fftshift)[freqstart:freqend] |> normalizer )
spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )


F = Figure(resolution=(800,600))
ax = Axis(F[1,1],
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(7),
yscale = log10)

#signal
scatter!(ax,freq_domain,spec_signal .*10^8,linewidth=lw,
color=col_signal,marker=:x,markersize=ms)
scatter!(ax,freq_domain,spec_signal .*10^12,linewidth=lw,
color=col_signal,marker=:x,markersize=ms,label="signal")

#ssa
lines!(ax,freq_domain,spec_ssa,linewidth=lw,
color=col_ssa,linestyle=:solid)
lines!(ax,freq_domain,spec_res_ssa,linewidth=lw,
color=col_ssa,linestyle=:solid)

#nlsa
lines!(ax,freq_domain,spec_nlsa,linewidth=lw,
color=col_nlsa,linestyle=:solid)
lines!(ax,freq_domain,spec_res_nlsa,linewidth=lw,
color=col_nlsa,linestyle=:solid)

hideydecorations!(ax)
hidespines!(ax, :t, :r, :l) # only top and right
#axislegend(ax)

save(dir*"comp_spec.png",F)















freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

Four = fft(signal) |> fftshift
signal_domain = abs.(Four[freqstart:freqend]) .* 1000 .+ 10^-6
Four = fft(trend) |> fftshift
trend_domain = abs.(Four[freqstart:freqend]) .+ 10^-6
Four = fft(signal .- trend) |> fftshift
res_domain = abs.(Four[freqstart:freqend]) ./ 1000 .+ 10^-6

ax_nlsa = Axis(F[1,1:2],title="the first yearly double mode captures $val",
yticksvisible = false,yticklabelsvisible = false,xlabel="t (d)",ylabel="GPP")
ax_nlsa_f = Axis(F[1,3],title="FFT spec",yscale=log10,
yticksvisible = false,yticklabelsvisible = false,xlabel="f (/a)",ylabel="log rel power")

lines!(ax_nlsa,1:N,signal .- trend .-4 ,color="red",label="residual")
lines!(ax_nlsa,1:N,signal,color="grey",label="signal")
lines!(ax_nlsa,1:N,trend,color="blue",label="trend")
axislegend(ax_nlsa)

lines!(ax_nlsa_f,freq_domain,res_domain,color="red")
lines!(ax_nlsa_f,freq_domain,signal_domain,color="grey")
lines!(ax_nlsa_f,freq_domain,trend_domain,color="blue")

scatter!(ax_nlsa_f,freq_domain,res_domain,color="red",marker='x')
scatter!(ax_nlsa_f,freq_domain,signal_domain,color="grey",marker='+')
scatter!(ax_nlsa_f,freq_domain,trend_domain,color="blue",marker='o')


Label(F[0, :], text = "comparison single mode and harmonics | $(meta[spot,"name"]) $(meta[spot,"IGBP_class"]) $(variable_list[vari]) $(method)", fontsize = 30)




freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

Four = fft(signal) |> fftshift
signal_domain = abs.(Four[freqstart:freqend]) .* 1000 .+ 10^-6
Four = fft(trend) |> fftshift
trend_domain = abs.(Four[freqstart:freqend]) .+ 10^-6
Four = fft(signal .- trend) |> fftshift
res_domain = abs.(Four[freqstart:freqend]) ./ 1000 .+ 10^-6

ax_ssa = Axis(F[2,1:2],title="$le identified trend modes capture $val | $val_m lost by mixing of $le_m modes",
yticksvisible = false,yticklabelsvisible = false,xlabel="t (d)",ylabel="GPP")
ax_ssa_f = Axis(F[2,3],title="FFT spec",yscale=log10,
yticksvisible = false,yticklabelsvisible = false,xlabel="f (/a)",ylabel="log rel power")

lines!(ax_ssa,1:N,signal .- trend .- 4,color="red")
lines!(ax_ssa,1:N,signal,color="grey")

lines!(ax_ssa,1:N,trend,color="blue")


lines!(ax_ssa_f,freq_domain,res_domain,color="red")
lines!(ax_ssa_f,freq_domain,signal_domain,color="grey")
lines!(ax_ssa_f,freq_domain,trend_domain,color="blue")

scatter!(ax_ssa_f,freq_domain,res_domain,color="red",marker='x')
scatter!(ax_ssa_f,freq_domain,signal_domain,color="grey",marker='+')
scatter!(ax_ssa_f,freq_domain,trend_domain,color="blue",marker='o')

save(savedir*"$(spot)_$(vari)_$(method).png",F)
#save(dir*"test.png",F)

end
