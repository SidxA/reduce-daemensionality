window_range = Array(120:120:3600)
k=32
p = plot(title="V")
for W in window_range
    pldiff = []
    plssa = []
    diff = npzread(dir*"runs/diffg_$W.npz")
    ssa = npzread(dir*"runs/ssag_$W.npz")
    diffHT = hcat([atan.(imag(hilbert_transform(diff[:,i])),real(hilbert_transform(diff[:,i]))) for i in 1:k]...)
    ssaHT = hcat([atan.(imag(hilbert_transform(ssa[:,i])),real(hilbert_transform(ssa[:,i]))) for i in 1:k]...)
    diffT = 19.55 ./ [count_sign_flip(diffHT[:,i]) for i in 1:k]
    ssaT = 19.55 ./ [count_sign_flip(ssaHT[:,i]) for i in 1:k]
    ssaV = vcat([var(ssa[:,i]) for i in 1:k]...)
    diffV = vcat([var(diff[:,i]) for i in 1:k]...)
    diff_s = 10^-6
    diff_ns = 10^-6
    ssa_s = 10^-6
    ssa_ns = 10^-6
    for i in 1:20
        if 0.45 < diffT[i] < 0.55
            diff_s += diffV[i]
        elseif diffT[i] > 0.55
            diff_ns += diffV[i]
        end

        if 0.45 < ssaT[i] < 0.55
            ssa_s += ssaV[i]
        elseif ssaT[i] > 0.55
            ssa_ns += ssaV[i]
        end
    end
    println(diff_ns)
    p = scatter!((W,diff_s),c = :red,marker=:X,markersize=3)
    p = scatter!((W,diff_ns),c = :red,marker=:o,markersize=3)
    p = scatter!((W,ssa_s),c = :blue,marker=:X,markersize=3)
    p = scatter!((W,ssa_ns),c = :blue,marker=:o,markersize=3)
end



for j in 1:20
    if diffT[j] >= 0.6
        p = scatter!((W,diffT[j]*2),c = :red,marker=:X,markersize=3)
    end
    if ssaT[j] >= 0.6
        p = scatter!((W,ssaT[j]*2),c = :blue,marker=:X,markersize=3)
    end
end

    #p = scatter(ssaT,ssaV,label="ssa",legend=:none,markersize=3)
    #p = scatter!(diffT,diffV,label="diff",xlabel="T",ylabel="var",markersize=3)
    #p = scatter!((W / 365,0.1),label="W",c=:black,marker=:X,legend=:none,xlim=(0,10),xaxis=(Array(1:1:10)))
    #pl = push!(pl,p)
    #savefig(plot(p),dir*"runs/hist_$W")

#=
@everywhere begin
    window_range = Array(120:120:3900)
    include("/net/home/lschulz/scripts/toolbox.jl")
    data = centralizer(load_data())

    function run_ssa(window_length)
        npzwrite(dir*"runs/ssag_$window_length.npz",SSA_g(data,window_length,32))
        return nothing
    end

    function run_diff(window_length)
        npzwrite(dir*"runs/diffg_$window_length.npz",Diff_g(data,window_length,32))
        return nothing
    end

end


seasons_gssa = sum(SSA_g(signal,window,20)[:,1:2],dims=2)
seasons_gdiff = sum(Diff_g(signal,window,20)[:,1:2],dims=2)


seasons_lssa1 = sum(local_SSA(signal,window,win2)[:,1:2],dims=2)
seasons_lssa2 = sum(local_SSA(signal,window,win2)[:,1:2],dims=2)
seasons_lssa3 = sum(local_SSA(signal,window,win2)[:,1:2],dims=2)

for (i,x) in enumerate([seasons_gdiff[:,1],seasons_gssa[:,1]])
    data = signal .- x
    ssa = SSA_g(data,3000,20)
    diff = Diff_g(data,3000,20)
    diffHT = hcat([atan.(imag(hilbert_transform(diff[:,i])),real(hilbert_transform(diff[:,i]))) for i in 1:20]...)
    ssaHT = hcat([atan.(imag(hilbert_transform(ssa[:,i])),real(hilbert_transform(ssa[:,i]))) for i in 1:20]...)
    diffT = 19.55 ./ [count_sign_flip(diffHT[:,i]) for i in 1:20]
    ssaT = 19.55 ./ [count_sign_flip(ssaHT[:,i]) for i in 1:20]
    ssaV = vcat([var(ssa[:,i]) for i in 1:20]...)
    diffV = vcat([var(diff[:,i]) for i in 1:20]...)
    p = scatter(ssaT,ssaV,label="ssa")
    p = scatter!(diffT,diffV,label="diff",xlabel="T",ylabel="var")
    savefig(plot(p),dir*"hist_$i")
    println(i)
end
=#