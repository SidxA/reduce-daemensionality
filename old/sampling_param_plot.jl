data = Array(CSV.read(dir*"sample_ssa_log.txt", DataFrame))[:,1:end-1]
samples = size(data)[1]

B_f = []
W_f = []
k_f = []
Ws_f = []
ks_f = []
B_h = []
W_h = []
k_h = []
Ws_h = []
ks_h = []

for i in 1:samples
    name = data[i,1]
    N = data[i,2]
    B = data[i,3]
    W = data[i,4]
    k = data[i,5]
    Ws = data[i,6]
    ks = data[i,7]
    Tf = filter(!iszero, data[i,8:8+k_max])
    Th = filter(!iszero, data[i,9+k_max:end])

    for T in Tf
        B_f = push!(B_f,(B,T))
        W_f = push!(W_f,(W,T))
        k_f = push!(k_f,(k,T))
        Ws_f = push!(Ws_f,(Ws,T))
        ks_f = push!(ks_f,(ks,T))
    end

    for T in Th
        B_h = push!(B_h,(B,T))
        W_h = push!(W_h,(W,T))
        k_h = push!(k_h,(k,T))
        Ws_h = push!(Ws_h,(Ws,T))
        ks_h = push!(ks_h,(ks,T))
    end
end

B_f =   [B_f[i] for i in 1:length(B_f)]
W_f =   [W_f[i] for i in 1:length(W_f)]
k_f =   [k_f[i] for i in 1:length(k_f)]
Ws_f =  [Ws_f[i] for i in 1:length(Ws_f)]
ks_f =  [ks_f[i] for i in 1:length(ks_f)]
B_h =   [B_h[i] for i in 1:length(B_h)]
W_h =   [W_h[i] for i in 1:length(W_h)]
k_h =   [k_h[i] for i in 1:length(k_h)]
Ws_h =  [Ws_h[i] for i in 1:length(Ws_h)]
ks_h =  [ks_h[i] for i in 1:length(ks_h)]

function plotsave_bothlist(list1,list2,name,)
    p = scatter(list1,c=:black,marker = :x,xlabel="parameter",ylabel="detected T")
    p = scatter!(list2,c=:red,marker = :x)
    savefig(plot(p,dpi=800,legend=:none),dir*name)
end

plotsave_bothlist(B_f,B_h,"B")
plotsave_bothlist(W_f,W_h,"W")
plotsave_bothlist(k_f,k_h,"k")
plotsave_bothlist(Ws_f,Ws_h,"Ws")
plotsave_bothlist(ks_f,ks_h,"ks")