
include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie
"""
yearly trends
"""


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
Ts = 1 / 365.25
# Start time 
t0 = 0 
tmax = t0 + (N-1) * Ts
# time coordinate
t = t0:Ts:tmax


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

signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]

matrix = hcat([signal[i] for i in year_ind]...)


f,ax,s = series(1:365,matrix',color=:tokyo,linewidth=1)

ssa_rc = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["RC"]
nlsa_rc = load(create_file_list(outdir,"diff",W,vari,preproc)[spot])["RC"]

ssa = sum(ssa_rc,dims=2)[:]
nlsa = sum(nlsa_rc,dims=2)[:]

ssa_traj = hcat([ssa[i] for i in year_ind]...)
nlsa_traj = hcat([nlsa[i] for i in year_ind]...)

gradi = :viridis
F = Figure()
ax_ssa = Axis(F[1,1])
ax_nlsa = Axis(F[1,2])

series!(ax_ssa,1:365,ssa_traj',color=gradi,linewidth=1)
series!(ax_nlsa,1:365,nlsa_traj',color=gradi,linewidth=1)

save(dir*"test.png",F)


"""
now the differences
"""

#savedirname = dir*"poster/"
spot = 3
threshold = 10
vari=1


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


colors = cgrad(:viridis,15,categorical=true)
years = 1985 .+ (1:15)
for spot = 1:20

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

# nlsa reconstructi

colors = cgrad(:viridis,15,categorical=true)
years = 1985 .+ (1:15)
for spot = 1:20

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


ssa_fund = hcat([ssa_trend_first[i] for i in year_ind]...)
ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
nlsa_fund = hcat([nlsa_trend_first[i] for i in year_ind]...)
nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)


f = Figure()
Axis(f[1, 1], title =  "in $spot",
    yticks = ((1:15) ./ 2 .-1,  string.(years)),
    xticks = ([365/2, 400 + 365/2, 800 + 365/2, 1200 + 465/2],
    ["SSA fund","SSA harm","NLSA fund","NLSA harm"])
)
for i in 1:15
    d = lines!(1:365,ssa_fund[:,i].+(0.5*i),color = colors[i])#,offset = i / 4,)
    d = lines!(400 .+ (1:365),ssa_harm[:,i].+(0.5*i),color = colors[i])
    d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
    d = lines!(1200 .+ (1:365),nlsa_harm[:,i].+(0.5*i),color = colors[i])
end

save(dir*"years/$spot.png",f)

end
