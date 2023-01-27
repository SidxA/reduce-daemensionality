using Distributed
include("/net/home/lschulz/scripts/toolbox.jl")
addprocs(8)
@everywhere include("/net/home/lschulz/scripts/toolbox.jl")

function data_name_directory(datas,names)
    li = []
    for index in 1:length(names)
        if isdir("/net/home/lschulz/logs/KW_25/runs/"*names[index]) ==false
            mkdir("/net/home/lschulz/logs/KW_25/runs/"*names[index])
        end
        li = push!(li,(datas[index],"/net/home/lschulz/logs/KW_25/runs/"*names[index]))
    end
return li
end

@everywhere function g_ssa_data_directory_W(data,directory,W)
    npzwrite(string(directory,"/ssa_g_",W),SSA_g(data,W)[4])        #check k
    return nothing
end

W_list = Int.(100:100:800)

li = data_name_directory(
    vcat([npzread("/net/home/lschulz/logs/KW_25/simple")[i,:] for i in 1:4],
    [npzread("/net/home/lschulz/logs/KW_25/overlay")[i,:] for i in 1:4],
    [npzread("/net/home/lschulz/logs/KW_25/noise")[i,:] for i in 1:4]),
    vcat([[name*string(i)] for i=1:4,name = ["simple","overlay","noise"]]...)
)

for j in li
    pmap(W -> g_ssa_data_directory_W(centralizer(j[1]),j[2],W), W_list)
end
rmprocs(8)