#read out the log_RC_runs
function readoff_all()
    li = []

    for W=Array(100:100:6400), B = 3000, k= 1:20
        li = push!(li,readout(npzread(string("/net/home/lschulz/logs/KW_20/runs/SSA-g_RC_W_",W,"_B_2000"))[:,k]))
    end
return hcat(li...)
end

function readoff_byK()
    li = zeros(20,3,64)
    for (i,W)=enumerate(Array(100:100:6400)), B = 3000
        for k = 1:20
            li[k,:,i] = readout(npzread(string("/net/home/lschulz/logs/KW_20/runs/diff_g_RC_W_",W,"_B_3000"))[:,k])
        end
    end
return li
end

#component investigation
function readout(component::Vector)
    std_c = std(component)

    flips = count_sign_flip(component)
    maxima= count_maxima(component)

    #phase = phases(analyticsignal(component))
    #flips_p = count_sign_flip(phase)
    #maxima_p = count_maxima(phase)
    return [std_c,flips,maxima]#,flips_p,maxima_p]
end

#rc = readoff_all()
#=
x = vcat([fill(i,20) for i in Array(200:200:4000)]...)
savefig(heatmap(rc),dir*"5c")
std_c = rc[:,1]
T_flip1 = rc[:,2]
T_max1 = rc[:,3]
T_flip2 = rc[:,4]
T_max2 = rc[:,5]

savefig(scatter(x,std_c,legend=:none,alpha=alpha,yaxis=:log10),dir*"std_c")
savefig(scatter(x,T_flip1,legend=:none,alpha=alpha,yaxis=:log10),dir*"T_flip1")
savefig(scatter(x,T_max1,legend=:none,alpha=alpha,yaxis=:log10),dir*"T_max1")
savefig(scatter(x,T_flip2,legend=:none,alpha=alpha,yaxis=:log10),dir*"T_flip2")
savefig(scatter(x,T_max2,legend=:none,alpha=alpha,yaxis=:log10),dir*"T_max2")

markersizes = std_c
alpha = 0.3

savefig(scatter(T_flip1,std_c,legend=:none,alpha=alpha,xlabel="T",ylabel="std",xlim=(0,400),ylim=(0,400)),dir*"T_flip1_x_std")

=#

function plotthem(ssa_l,ssa_g,diff_g)
    a = 0.1
    m = 2
    s = :o

    p1=plot(xaxis=:log10,title="flip spectogram",legend=:none,xlabel="days period",ylabel="std",dpi=800)
    p2=plot(xaxis=:log10,title="maxima spectogram",legend=:none,xlabel="days period",ylabel="std",dpi=800)

    p1=scatter!(7304 ./ssa_l[:,2,:],ssa_l[:,1,:],marker=s,c=:green,alpha=a,markersize=m)
    p1=scatter!(7304 ./ssa_g[:,2,:],ssa_g[:,1,:],marker=s,c=:red,alpha=a,markersize=m)
    p1=scatter!(7304 ./diff_g[:,2,:],diff_g[:,1,:],marker=s,c=:blue,alpha=a,markersize=m)

    p2=scatter!(7304 ./ssa_l[:,3,:],ssa_l[:,1,:],marker=s,c=:green,alpha=a,markersize=m)
    p2=scatter!(7304 ./ssa_g[:,3,:],ssa_g[:,1,:],marker=s,c=:red,alpha=a,markersize=m)
    p2=scatter!(7304 ./diff_g[:,3,:],diff_g[:,1,:],marker=s,c=:blue,alpha=a,markersize=m)

    savefig(p1,dir*"results_f")
    savefig(p2,dir*"results_m")
    return nothing
end