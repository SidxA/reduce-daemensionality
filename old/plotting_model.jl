#=
function parameters()
    parameterset = []
    i=1
    for epsilon = [0.0,0.1,0.5]
        for W = [20,40,100,200,400]

            p = 0
            w = 2.5 * 10^-1
            parameterset = push!(parameterset,[i,p,epsilon,w,W])
            i+=1

            p = 1
            w = 2.5 * 10^-5
            parameterset = push!(parameterset,[i,p,epsilon,w,W])
            i+=1
            p = 2
            w = 2.5 * 10^-9
            parameterset = push!(parameterset,[i,p,epsilon,w,W])
            i+=1
        end
    end
    println("total parameter combinations $i")
    return parameterset
end


phases          = vcat([npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_phases") for i in 1:45]...)
protophases     = vcat([npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_protophases") for i in 1:45]...)
variances       = vcat([npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_variances") for i in 1:45]...)

#there are 45 parameters * 4 scenarios * 7 signals
scen_list = [plot(),plot(),plot(),plot()]
for noise_index = 1:3
    for scen_index = 1:4
            paindex = noise_index*.Array(1:15)
            phase = phases[paindex*scenindex*1,:,:]
            gSSA = cumsum(vcat([paindex*scenindex*protophases[(2:3)*k,:,:].+pi for k in 1:4]...),dims = 3)
            gDiff = cumsum(vcat([paindex*scenindex*protophases[(5:6)*k,:,:].+pi for k in 1:4]...),dims = 3)
            
            scenlist[scen_index] = plot!(plotlist,layout=(3,2))

parameterlist = parameters()
for para in parameterlist
    for scenario = 1:4    
    i = Int(para[1])
    p = para[2]
    epsilon = para[3]
    w = para[4]
    W = Int.(para[5])

    phases          = npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_phases")
    protophases     = npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_protophases")
    variances       = npzread("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_variances")

    phase = phases[1,:,:]


    gSSA = cumsum(vcat([protophases[(2:3)*k,:,:].+pi for k in 1:4]...),dims = 3)
    gDiff = vcat([protophases[(5:6)*k,:,:].+pi for k in 1:4]...)

    SSAplot = scatter(phase',legend=:none,markersize=3,c=:black)
    #SSAplot = scatter!(gSSA',legend=:none,markersize=1,c=:blue)

    Diffplot = scatter(phase',legend=:none,markersize=3,c=:black)
    #Diffplot = scatter!(gSSA',legend=:none,markersize=1,c=:blue)




#individual_plots = [plot(protophases[i,1,:],
#                    legend=:none,#string.(round.(variances[i,:]',digits=3)),
#                    title="$i",xticks=:none,yticks=:none,ylim=(0,7),xlim=(0,5000)) for i=1:4*(2*k+1)]
#savefig(plot(individual_plots...,layout=(4,2*k+1),dpi=1000),dir*"protophases")

vcat(npzread("/net/home/lschulz/logs/KW_25/runs/model/1_protophases"),
npzread("/net/home/lschulz/logs/KW_25/runs/model/2_protophases"))

=#

scenario_list = ["plain","centered","40 day cutoff","40 day cutoff + center"]
nsce    = 4
nnoise  = 3
np      = 3
nW      = 5
nN      = 5
nC      = 3
nT      = 8000
r_phases = Array{Float64}(undef,nsce,nnoise,np,nW,nN,nT)
r_SSA = Array{Float64}(undef,nsce,nnoise,np,nW,nN,nC,nT)
r_Diff = Array{Float64}(undef,nsce,nnoise,np,nW,nN,nC,nT)

#original data size = Nsce * Nnoi * Np * NW  * (1+Nco) x nN x nT
    for scen_index = 1:nsce
        for noise_index = 1:nnoise
            index = scen_index*noise_index*.*Array(1:15)


            phases_l = hcat([phases[i,:] for i in index]...)
            gSSA = cumsum(vcat([protophases[i*(2:3),:,:].+pi for i in index]...),dims = 3)
            gDiff = cumsum(vcat([protophases[i*(5:6),:,:].+pi for i in index]...),dims = 3)


            for n = 1:1
                r_phases[scen_index,noise_index,p_index,W_index,N_index,C_index,:]
                r_SSA[scen_index,noise_index,p_index,W_index,N_index,C_index,:]
                r_Diff[scen_index,noise_index,p_index,W_index,N_index,C_index,:]

            end




            
        ssaplot = scatter(phases_l,legend=:none,marker=:x,c=:black,markersize=1)
        diffplot = scatter(phases_l,legend=:none,marker=:x,c=:black,markersize=1)

            #one plot per scenario level
            p = plot(title="$(scenario_list[scen_index])")
            noise_list = [plot(),plot(),plot()]

            noise_list[noise_index] = plot([ssaplot,diffplot]...,layout=(1,2),title="noise $([0.0,0.1,0.5][noise_index])")
        p = plot!(noise_list...,layout=(3,1),dpi=800)
        savefig(plot(p),dir*"scenario_$scen_index")

        ssaplot = scatter!(gSSA[:,n,:]',legend=:none,marker=:x,c=:red,markersize=1)
        diffplot = scatter!(gDiff[:,n,:]',legend=:none,marker=:x,c=:blue,markersize=1)