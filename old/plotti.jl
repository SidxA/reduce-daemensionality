using FileIO
using ImageIO

img = FileIO.load("logs/KW_29/"*"worldmap.png")
table = CSV.read("logs/KW_29/"*"fluxsites.csv",DataFrame,header=0)
nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")

need to put yflip = false !

#lattitude is 6
# - longitude is 5

i = [findfirst(x->x==name,table[:,1]) for name in nameslist[1:end-1]]
i = [isnothing(ix) ? 0 : ix for ix in i]
i = i[findall(!iszero,i)]

p = plot()
p = plot!([-180,180],[90,-90],img,yflip=true)
p = scatter!(table[:,6],-table[:,5],alpha=0.2)
p = scatter!([(table[i,6],-table[ix,5]) for ix in i])
savefig(p,dir*"map")

#with europe
p = plot(legend=false)
p = plot!([-180,180],[90,-90],img,yflip=true)
p = scatter!([(table[i,6],-table[ix,5]) for ix in i],xlim=(-60,40),ylim=(-65,20))
savefig(p,dir*"europemap")

#between nodes
line(p,ind1,ind2,alpha) = plot!([table[i,:][ind1,6],table[i,:][ind2,6]],[-table[i,:][ind1,5],-table[i,:][ind2,5]],c=:black,linewidth=3,alpha=alpha)


"""
mapping to a picture of europe with slightly tweaked locations
"""
function lat_long()#return x,y
    #first use original locations
    table = CSV.read("logs/KW_29/"*"fluxsites.csv",DataFrame,header=0)
    filelist = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")[5:end-3]
    nameslist = String[]
    for file in filelist
        filenamesplit = split(file,"_")
        years = split(filenamesplit[end-1],"-")
        period = parse(Int,years[2])-parse(Int,years[1])
        if filenamesplit[4] == "ERAI" && filenamesplit[5] == "DD" && filenamesplit[end] =="beta-3.csv" && period > 30
        nameslist = push!(nameslist,filenamesplit[2])
        end
    end

    loc_ind=1:71
    i = [findfirst(x->x[1:2]==name[1:2],table[:,1]) for name in nameslist]
    i = [isnothing(ix) ? 0 : ix for ix in i]
    #israel would be around 15,40
    #island would be around 5, 50
    i[55] = 121
    i[56] = 27
    #loc_ind = loc_ind[findall(!iszero,i)]
    i = i[findall(!iszero,i)]

    #the stupid manual list
    tsml = 
    reshape(
        vcat(
    [
        #belgium
        [0,0+2],
        [2,0+2],
        [0,2+2],
        [-2,0+2],
        [0,-2+2],
        [-2,2+2],
        #switz?
        [0+2,0-2],
        [2+2,0-2],
        [0+2,2-2],
        [-2+2,0-2],
        [0+2,-2-2],
        [2+2,2-2],
        #czesk
        [0+2,0-2],
        [2+2,0-2],
        [0+2,2-2],
        [-2+2,0-2],
        [0+2,-2-2],
        [-2+2,2-2],
        #ger
        [0-2,0-1],
        [2-2,0-1],
        [0-2,2-1],
        [-2-2,0-1],
        [0-2,-2-1],
        [-2-2,2-1],
        [4-2,2-1],
        [4-2,0-1],
        [0-2,-4-1],
        [-4-2,0-1],
        [-4-2,-4-1],
        [-2-2,-4-1],
        #denmark
        [-2,-2+7],
        [-3,-4+7],
        #spain
        [-1,3],
        [-5,3],
        [-3,5],
        [-5,3],
        [-3,2],
        [-5,5],
        #finland
        [0,0],
        [2,0],
        [0,2],
        [-2,0],
        [0,-2],
        [-2,-2],
        #france
        [0,0-3],
        [2,0-3],
        [0,2-3],
        [-2,0-3],
        [0,-2-3],
        [-2,2-3],
        [2,2-3],
        [3,0-3],
        #gf
        [+45,+30],
        #gl
        [+45,-6],
        #ie
        [+10,0],
        #il
        [-5,+15],
        #italy
        [0,0],
        [2,0],
        [0,2],
        [-2,0],
        [0,-2],
        [-2,2],
        [4,2],
        [4,0],
        #russland
        [-2,-2],
        [-3,-4],
        #SE
        [0-4,0-2],
        [2-4,0-2],
        [0-4,2-2],
        [-2-4,0-2],
        [0-4,-2-2],
    ]...),
    2,71)'

    x = table[i,6]
    y = table[i,5]
    x = x .+ tsml[:,1]
    y = y .+ tsml[:,2]
    return x,y
end
"""
indices pairs of links that are inside k and epsilon
"""
using NearestNeighbors

function find_links(x,y,k,r)
    #neigbors
    data = hcat(x,y)'
    tree = KDTree(data)
    idxs_k,dists = knn(tree,data,k+1,true)
    #OR for balls of radius r
    idxs_r = inrange(tree, data, r)
    L = []
    for loc_ind in 1:71
        arr = intersect(idxs_k[loc_ind],idxs_r[loc_ind])
        self = findall(x->x==loc_ind,arr)
        arr = deleteat!(arr,self)
        for x in arr
            if isempty(findall(y->y==[x,loc_ind],L)) 
                L = push!(L,[loc_ind,x])
            end
        end
    end
    return Matrix(hcat(L...)')
end

"""
execution
"""
using CairoMakie
using MathTeXEngine
using FileIO
using ImageIO
img = FileIO.load("logs/KW_29/"*"worldmap.png")

x,y =lat_long()
x .+= rand(71).-0.5
y .+= rand(71).-0.5




r = 15
k = 7
F = Figure()
a = Axis(F[1,1],limits=(-10, 30, 30, 70))
image!(a,[-180,180],[-90,90],rotr90(img),yflip=true)
links = find_links(x,y,k,r)
start = links[:,1]
stop = links[:,2]
for l in 1:size(links)[1]
    lines!(a,[x[start][l],x[stop][l]],[y[start][l],y[stop][l]],color=:black)
end

scatter!(a,vec(x),vec(y))

save(dir*"map.png",F)

