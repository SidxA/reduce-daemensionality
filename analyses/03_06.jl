"""
recreate the bigf figure
"""


spot = 19
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