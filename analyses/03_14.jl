
using Pkg
Pkg.activate("Fouriertest")

using Dates
function init_logging()
    weekno = week(unix2datetime(time()))
    datestring = string("KW_2_",lpad(weekno,2,"0"),"/")
    workdir = "/net/home/lschulz/logs/"
    dir = workdir*datestring
    if isdir(dir)==false
        mkdir(dir)
    end
    return dir
end
dir = init_logging()


using Fouriertools
using DataFrames
using JLD2
using CairoMakie
using LsqFit

"""
all the plots
"""

#uses these functions

function reconstructor(pc,eof,N,W)
    P = N-W+1
    #A is the k_th PC (projection) of length P
    #rho is the k_th eigenvector (EOF) of length M

    #M is window length W
    R(t,M_t,L_t,U_t) = M_t*sum([pc[t-j+1]*eof[j] for j in L_t:U_t])
    function choicer(t)
        if 1<=t<=W-1
            return 1/t,1,t
        elseif W<=t<=P
            return 1/W,1,W
        elseif P+1<=t<=N
            return 1/(N-t+1),t-N+W,W
        end
    end
    return [R(t,choicer(t)...) for t in 1:N]
end


function create_file_list(outdir,method,W,vari,preproc) #at the moment this does diffusion
    allfiles = readdir(outdir)
    filenamesplit = split.(allfiles,"_")
    method_inds = findall(x->( x[1]==method && x[2] == "$W" && x[4] == "$vari" && x[end][1:4] == preproc),[filenamesplit[i] for i in 1:length(filenamesplit)])
    files = outdir .* allfiles[method_inds]
    return files
end #return files


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

#prerequisites
#spot = 1
outdir="/net/scratch/lschulz/fluxfullset_1a/"
#meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 365
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
tw = t0:Ts:(t0+(W-1)*Ts)

color_ssa = "darkgreen"
color_nlsa = "purple"
color_signal = "grey50"
lw = 1
ms = 5

years = ((1:5478) ./ 365) .+ 2005

# plotting parameters

meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]

igbp = meta[spot,"IGBP_class"]
name = meta[spot,"name"]
variable = ["GPP"][vari]

titlestring = "$variable in $spot $name $igbp"
savedirstring = "$(variable)_$(spot)_$(W)d_"
savedir = savedirname*savedirstring

#load files
Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
file_ssa = load(Filename_ssa)
rc_ssa = sum(file_ssa["RC"],dims=2)[:]
file_nlsa = load(Filename_nlsa)
rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]

signal = file_ssa["signal"]

#FFT stuff
ms = 8
odi = 4.5 #log offset 
normalizer(x) = x./maximum(x)
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]


spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )
spec_ssa = (abs.(fft(rc_ssa) |> fftshift)[freqstart:freqend] |> normalizer )
spec_nlsa = (abs.(fft(rc_nlsa) |> fftshift)[freqstart:freqend] |> normalizer )

"""
build the 1 year spectrum peak
"""
    Eof = load(Filename_ssa)["EOF"]
    threshold = 10

    W,k = size(Eof)

    for i in 1:k

        mode = Eof[:,i]

        f = Figure()
        ax = Axis(f[1,1],limits=(0,4,0,20),xticks=0:0.5:4,)

        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1/Ts) |> fftshift




        #does resize it, but spectrum should be intact !

        mode2 = resample(mode,W*2,normalize=true)
        

        Four2 = abs.(fft(mode2) |> fftshift)
        freqs2 = fftfreq(2*length(tw), 2/Ts) |> fftshift

        Four2 ./= maximum(Four2)

        #need to create a mask to isolate the data between 0 and 8

        f0 = findall(x->(abs(x-0)<=0.01),freqs2)[end]
        f8 = findall(x->(abs(x-8)<=0.01),freqs2)[1]

        Four2[1:f0] .= 0
        Four2[f8:end] .= 0

        f,v,s = fit_gauss(freqs2,abs.(Four2))
        peak = gauss(freqs2,[f,v,s])

        F = Figure()
        ax = Axis(F[1,1],limits=(0,8,0,1),xticks=0:0.5:8)

        scatter!(ax,freqs2,abs.(Four2))
        scatter!(ax,freqs2,peak)

        freqstart2 = findall(x->x>=1/12,freqs2)[1]
        freqend2 = findall(x->x>=6,freqs2)[1]

        
        val_domain = abs.(Four2[freqstart:freqend])
        #val_domain ./= maximum(val_domain)

        freq_domain = freqs2[freqstart:freqend]





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


# ssa reconstruction: freq harmonics

file = load(Filename_ssa)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]

li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
ssa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

# nlsa reconstruction: freq harmonics

file = load(Filename)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]

li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
nlsa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

#the fundamental

li_harmonics=[1,2]
val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
#val_m = (round(sum(lambda[li_mixed]),digits=3))
#le_m = length(li_mixed)
trend_first = sum(RC[:,li_harmonics],dims=2)[:]

#the figure

offset = 2.5
lw = 3
lw_s = 1
ms = 4

F = Figure(resolution=(800,800))


ax_time = Axis(F[1,1],
xticks = 2005:3:2020,
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(3),
xlabel = "time (a)",
ylabel = "GPP"
)


#signal background

scatter!(ax_time,years,signal .+ offset,linewidth=lw,
color=color_signal,marker=:x,markersize=ms)
lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
color=color_signal,marker=:x,markersize=ms)
#scatter!(ax_time,years,signal .+ 12.0,linewidth=lw,
#color=color_signal,marker=:x,markersize=ms,label="signal")

#ssa
lines!(ax_time,years,ssa_trend_harm .+ offset,linewidth=lw,
color=color_ssa,linestyle=:solid,label="ssa")
lines!(ax_time,years,signal .- ssa_trend_harm .+ 0,linewidth=lw_s,
color=color_ssa,linestyle=:solid,label="ssa")

#nlsa
lines!(ax_time,years,nlsa_trend_harm .+ offset,linewidth=lw,
color=color_nlsa,linestyle=:solid,label="nlsa")
lines!(ax_time,years,signal .- nlsa_trend_harm,linewidth=lw_s,
color=color_nlsa,linestyle=:solid,label="nlsa")

#hideydecorations!(ax_time)
#hidespines!(ax_time, :t, :r, :l) # only top and right
#axislegend(ax_time)

#save(savedir*"$(spot)_comp.png",F)

#now the additional FFT

ms = 8
odi = 4.5
normalizer(x) = x./maximum(x)
freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

spec_ssa = (abs.(fft(ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^odi
spec_res_ssa = (abs.(fft(signal .- ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer)
spec_nlsa = (abs.(fft(nlsa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^odi
spec_res_nlsa = (abs.(fft(signal .- nlsa_trend_harm) |> fftshift)[freqstart:freqend] |> normalizer )
spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )


ax_spec = Axis(F[2,1],
xminorticksvisible = true,
xminorgridvisible = true,
xminorticks = IntervalsBetween(7),
yscale = log10,
xlabel = "frequency (/a)",
ylabel = "relative power")

#signal
#scatter!(ax_spec,freq_domain,spec_signal,linewidth=lw,
#color=color_signal,marker=:x,markersize=ms)
scatter!(ax_spec,freq_domain,spec_signal .*10^odi,linewidth=lw,
color=color_signal,marker=:x,markersize=ms,label="signal")
#lines!(ax_spec,freq_domain,spec_signal,linewidth=lw_s,
#color=color_signal,linestyle=:solid)
lines!(ax_spec,freq_domain,spec_signal .* 10^odi,linewidth=lw_s,
color=color_signal,linestyle=:solid)

#ssa
lines!(ax_spec,freq_domain,spec_ssa,linewidth=lw,
color=color_ssa,linestyle=:solid)
lines!(ax_spec,freq_domain,spec_res_ssa,linewidth=lw,
color=color_ssa,linestyle=:solid)

#nlsa
lines!(ax_spec,freq_domain,spec_nlsa,linewidth=lw,
color=color_nlsa,linestyle=:solid)
lines!(ax_spec,freq_domain,spec_res_nlsa,linewidth=lw,
color=color_nlsa,linestyle=:solid)

#hideydecorations!(ax_spec)
#hidespines!(a_spec, :t, :r, :l) # only top and right
#axislegend(ax_spec)

Label(F[0, 1], text = "seasons + FFT for "*titlestring,
fontsize = 24)

#trim!(F.layout)

save(savedir*"seasonal.png",F)

"""
the run on the years
"""


year_ind =
[
1:365,
366:730,
731:1095,
1096:1460,
1462:1826, #89 1 offset
1827:2191, 
2192:2556,
2557:2921,
2923:3287, # 93 1 offset
3288:3652,
3653:4017,
4018:4382,
4384:4748, #97 1 offset
4749:5113,
5114:5478,
]

years_plot = 2005 .+ (1:15)

color_signal_grad =  cgrad([:grey20,:grey80],15,categorical=true,rev=true)
color_ssa_grad =  cgrad([:darkgreen,:lightgreen],15,categorical=true,rev=true)
color_nlsa_grad =  cgrad([:purple,:plum],15,categorical=true,rev=true)

lw = 1


fund = hcat([trend_first[i] for i in year_ind]...)

ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)


f = Figure(resolution=(800,800),)

ax_years = Axis(f[1, 1],
    yticks = ((1:15) ./ 2 .-1,  string.(years_plot)),
    xticks = ([365/2, 400 + 365/2, 800 + 365/2],
    ["fund","SSA harm","NLSA harm"]),
    ylabel="years"
)


for i in 1:15
    d = lines!(ax_years,1:365,fund[:,i].+(0.5*i),color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
    d = lines!(ax_years,400 .+ (1:365),ssa_harm[:,i].+(0.5*i),color = color_ssa_grad[i],linewidth=lw)
    #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
    d = lines!(ax_years,800 .+ (1:365),nlsa_harm[:,i].+(0.5*i),color = color_nlsa_grad[i],linewidth=lw)
end

ax_years_2 = Axis(f[2,1],
xticks = ([365/2, 400 + 365/2, 800 + 365/2],
    ["fund","SSA harm","NLSA harm"]),
    ylabel="GPP"
    )

for i in 1:15
    d = lines!(ax_years_2,1:365,fund[:,i],color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
    d = lines!(ax_years_2,400 .+ (1:365),ssa_harm[:,i],color = color_ssa_grad[i],linewidth=lw)
    #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
    d = lines!(ax_years_2,800 .+ (1:365),nlsa_harm[:,i],color = color_nlsa_grad[i],linewidth=lw)
end

Label(f[0, :], text = "changes in seasonal behavior for "*titlestring,
fontsize = 24)

#trim!(f.layout)

save(savedir*"years.png",f)
