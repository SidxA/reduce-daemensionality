#=
include("/net/home/lschulz/scripts/toolbox.jl")
names = readdir(dir*"TA_ERA_DAY")
runs = dir*"TA_ERA_DAY_runs/"

name = names[rand(1:72)]

signal = centralizer(npzread(dir*"TA_ERA_DAY/"*name))

T = length(signal)
years = T / 365
k = 16

diff = Diff_g(signal,W,k)
diffHT = hcat([atan.(imag(hilbert_transform(diff[:,i])),real(hilbert_transform(diff[:,i]))) for i in 1:k]...)
diffT = [count_sign_flip(diffHT[:,i]) for i in 1:k]
diffV = vcat([var(diff[:,i]) for i in 1:k]...)

eps = 0
for i in 1:k
    t = diffT[i]
    v = diffV[i]
    if norm(t - T/365 *2) > eps && norm(t - T/365 *3) > eps && norm(t - T/365 *4) > eps && v > 10^-4
        println("Diff $t \t $(round(T/t *2 / 365.25 * 12,digits = 2)) m \t \t $v")
    end
    t = ssaT[i]
    v = ssaV[i]
    if norm(t - T/365 *2) > eps && norm(t - T/365 *3) > eps && norm(t - T/365 *4) > eps && v > 10^-4
        println("SSA $t \t $(round(T/t *2 / 365.25 * 12,digits = 2)) m \t \t $v")
    end
end

ssa = SSA_g(signal,W,k)
ssaHT = hcat([atan.(imag(hilbert_transform(ssa[:,i])),real(hilbert_transform(ssa[:,i]))) for i in 1:k]...)
ssaT = [count_sign_flip(ssaHT[:,i]) for i in 1:k]
ssaV = vcat([var(ssa[:,i]) for i in 1:k]...)

if norm(t - T/365 *2) > 1 || norm(t - T/365 *3) > 1 || norm(t - T/365 *2) > 4
    println("$t \t $v")
end

=#



#=
function  mode_extraction(signal,W,k)

    N=length(data)
    P = N - W +1
    emb_data = embed_lag(data,W)

    pca = fit(PCA,emb_data,method=key,pratio=1.0,maxoutdim=k)
    EOFssa = projection(pca)

    t=1
    α=0.5
    ɛ=512

    diff = fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=α, ɛ=ɛ)
    EOFdiff = ManifoldLearning.transform(diff)'





something like an automated plot procedure?
something about how much this periodicity captures ... like how much variance does the rest still have?
and how much variance is in the seasons


savefig(plot([diffi[:,4].+diffi[:,5].+0.2 
diffi[:,6].+diffi[:,7].-0.2 
diffi[:,8].+diffi[:,9].+0.4 
diffi[:,10].+diffi[:,11] 
diffi[:,12].+diffi[:,13].-0.4 
diffi[:,14].+diffi[:,15].-0.6],
label=[182 121 91 73 60 52]),dir*"diff")

=#
#epsi_T is samples difference in periodlength that still count as the same, eg. 2,3,4
#epsi_A is the required crosscorrelation at the maxima in order to count
#returns pairs, periods,solo

using Distributed

addprocs(75)
@everywhere begin

    include("/net/home/lschulz/scripts/toolbox.jl")
    Lnames = readdir("/net/home/lschulz/logs/KW_26/TA_ERA_DAY")

    
    #W = 4500
    k = 75
    W_L = Array(500:25:5800)
    function param()
        parada = []
        i=1
        for t in [1,2,3]
            for alpha in [0,0.25,0.5,0.75,1]
                for epsi in [64,128,256,512]
                    parada = push!(parada,Float64.([t,alpha,epsi,i]))
                    i+=1
                end
            end
        end
        return parada
    end

    function paint(para)
        t=Int.(para[1])
        α=para[2]
        ɛ=Int.(para[3])
        i = Int.(para[4])

        N=length(signal)
        P = N - W +1

        #mo,per,tre = modes(signal,ManifoldLearning.transform(fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=α, ɛ=ɛ))')
        npzwrite(dir*"diff/EOF_$i",ManifoldLearning.transform(fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=α, ɛ=ɛ))')
        println("$i")
    end

    function SSA_periods(W)
        data = centralizer(npzread("/net/home/lschulz/logs/KW_26/TA_ERA_DAY/"*Lnames[62]))
        N=length(data)
        #W = Int(floor(N/2))-5
        P = N - W +1
        emb_data = embed_lag(data,W)
        pca = fit(PCA,emb_data,method=:auto,pratio=1.0,maxoutdim=k)
        EOF = projection(pca)
        npzwrite(dir*"W/ssaFull_W_$W",EOF)
        #=
        PC = hcat([pc!(data,EOF[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOF[:,i],N,W) for i in 1:k]...)
        for p = [(1,0.7),(5,0.7),(10,0.7),(1,0.1),(5,0.1),(10,0.1)]
            epsi_T = p[1]
            epsi_A = p[2]
            pairs,periods,solos = pairing(EOF,epsi_T,epsi_A)
            n_pairs = length(pairs)
            n_solo = length(solos)
            modi = hcat([RC[:,pair[1]].+RC[:,pair[2]] for pair in pairs]...)
            trends = hcat([RC[:,solo] for solo in solos]...)
            if isempty(trends)
                trends = zeros(N)
            elseif isempty(modi)
                modi = zeros(N)
            end
            npzwrite(dir*"W/ssaFull_W_$(W)_$p",hcat(modi,trends))
            println("$W \t $p \t $(periods[periods .> 370.0])")
        end
        =#
        return nothing
    end

    function SSA_outlet(W)
        EOF = npzread(dir*"W/ssaFull_W_$W")
        epsi_T = 30
        epsi_A = 0.1
        pairs,periods = pairing_full(EOF,epsi_T,epsi_A,false)
        per_cc = periods[periods .> 500.0]
        k =size(EOF)[2]
        N = 11688
        pe = Float64[0.5 .* N ./ iterated_hilbert(EOF[:,i],64)[1] for i=1:k]
        per_hht = pe[pe .> 500.0]
        println("$W \t CC $(per_cc) \t HHT $(per_hht)")
        return nothing
    end

    function iterated_hilbert(mode,l)
        ll = Int64[]
        for i = 1:l
            HilbertT = hilbert_transform(mode)
            protophases = atan.(imag(HilbertT),real(HilbertT))
            s = count_sign_flip(protophases)
            #println(s)
            ll = append!(ll,s)
            mode = protophases
        end
        return median(ll),var(ll)/sqrt(l)
    end

    function pairing(EOFs,epsi_T,epsi_A,print)
        pairs = []
        solo = []
        periods = []
        k = size(EOFs)[2]
        indices = Array(1:k)
        lags = 1:Int(size(EOFs)[1]-1)
        i = 1
        j = 2
        while isempty(indices) == false
            if length(indices) == 1
                solo = push!(solo,indices[1])
                break
            end
            if print println("indices $indices") end
            x=EOFs[:,i]
            y=EOFs[:,j]
            c = crosscor(x,y,lags)
            cc = diff(argmaxima(c))
            if length(cc) < 3
                if print println("$i \t $j less then 3 correlation maxima") end
                if j == indices[end]
                    solo = push!(solo,i)
                    indices = deleteat!(indices,findall(x->x==i,indices))
                    i = indices[1]
                    j = indices[2]
                else
                    j = indices[findall(x->x==j,indices)[1]+1]
                end
            elseif length(indices) == 2
                if abs(cc[1]-cc[2])<=epsi_T && all(c[argmaxima(c)][1:2] .> epsi_A)
                    pairs = push!(pairs,(i,j))
                    periods = push!(periods,mean(cc[1:4]))
                    break
                else
                    solo = push!(solo,indices[1])
                    solo = push!(solo,indices[2])
                    break
                end
            elseif abs(cc[1]-cc[2])<=epsi_T && all(c[argmaxima(c)][1:2] .> epsi_A)
                if print println("$i \t $j found T $(cc[1])") end
                pairs = push!(pairs,(i,j))
                periods = push!(periods,mean(cc[1:2]))
                indices = deleteat!(indices,findall(x->x==i,indices))
                indices = deleteat!(indices,findall(x->x==j,indices))
                if length(indices) == 1
                    solo = push!(solo,indices[1])
                    break
                end
                i = indices[1]
                j = indices[2]
            elseif j == indices[end]
                solo = push!(solo,i)
                indices = deleteat!(indices,findall(x->x==i,indices))
                i = indices[1]
                j = indices[2]
            else
                if print println("$i \t $j not inside scope") end
                j = indices[findall(x->x==j,indices)[1]+1]
            end
        end
        return pairs, periods,solo
    end

    function pairing_full(EOFs,epsi_T,epsi_A,print)
        pairs = Tuple[]
        periods = Float64[]
        k = size(EOFs)[2]
        indices = Array(1:k)
        lags = 1:Int(size(EOFs)[1]-1)
        for i=1:k, j= 1:k
            x=EOFs[:,i]
            y=EOFs[:,j]
            c = crosscor(x,y,lags)
            cc = diff(argmaxima(c))
            if length(cc) >= 2
                if abs(cc[1]-cc[2])<=epsi_T && all(c[argmaxima(c)][1:2] .> epsi_A)
                    if print println("$i \t $j found T $(cc[1])") end
                    pairs = push!(pairs,(i,j))
                    periods = push!(periods,mean(cc[1:2]))
                end
            end
        end
        return pairs,periods
    end

    #returns modes,periods,trends
    function modes(signal,EOF)
        W = size(EOF)[1]
        k = size(EOF)[2]
        N = length(signal)
        P = N - W +1

        epsi_T = 3
        epsi_A = 0.5

        pairs,periods,solos = pairing(EOF,epsi_T,epsi_A)
        n_pairs = length(pairs)
        n_solo = length(solos)
        PC = hcat([pc!(signal,EOF[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOF[:,i],N,W) for i in 1:k]...)

        modes = hcat([RC[:,pair[1]].+RC[:,pair[2]] for pair in pairs]...)
        trends = hcat([RC[:,solo] for solo in solos]...)
        return modes,periods,trends
    end

end

pmap(p->SSA_periods(p),W_L)
rmprocs(75)
#=
#p_alpha = plot()
#p_epsilon = plot()
#p_t = plot()
#p = param()
#epsi_T = 4
epsi_A = 0.4
for i = 1:60
    t = p[i][1]
    a = p[i][2]
    e = p[i][3]
    pairs,periods,solo = pairing(l[i],epsi_T,epsi_A)
    if isempty(periods) == false
        #p_alpha = scatter!(a.*ones(length(periods)),periods)
        p_t = scatter!(t.*ones(length(periods)),periods)
        #p_epsilon = scatter!(e.*ones(length(periods)),periods)
    end

end
savefig(plot(p_t,legend=:none),dir*"t_t")
=#