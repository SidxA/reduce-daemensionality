"""

TA temperature air,
PA pressure air,
VPD vapur pressure deficit,
WS wind speed/ WD direction,
P precipitation,
SW_IN shortwave,
LW_IN longwave,
LW_IN_JSB

NETRAD 		Net radiation
PPFD_IN 		Photosynthetic photon flux density, incoming
PPFD_OUT 		Photosynthetic photon flux density, outgoing

NET ECOSYSTEM EXCHANGE
respiration
GPP


TA_F_MDS_DAY 		Average daytime TA_F_MDS
DD 	deg C 	from half-hourly data

TA_F_MDS_DAY_SD 		Standard deviation for TA_F_MDS_DAY
DD 	nondimensional 	fraction between 0-1, indicating percentage of measured and good quality gapfill data

TA_ERA_DAY 		Average daytime TA_ERA


TA_F 		Air temperature, consolidated from TA_F_MDS and TA_ERA


TA_F_QC 		Quality flag for TA_F
DD 	nondimensional 	fraction between 0-1, indicating percentage of measured and good quality gapfill data

TA_F_DAY 		Average daytime TA_F

ls /net/data/LAI/LAI_AVHRR contains almost 40 years Leaf Area Index

"""

# choose erai fluxnet 

filelist = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")[5:end-3]
put_in = "/net/scratch/lschulz/FLUXERAI/"

#checked that tim eis correct!
N = 11322
i=1
j=N

TIMESTAMP = []
TA_ERA = []
TA_ERA_NIGHT = []
TA_ERA_NIGHT_SD = []
TA_ERA_DAY = []
TA_ERA_DAY_SD = []
SW_IN_ERA = []
LW_IN_ERA = []
VPD_ERA = []
PA_ERA = []
P_ERA = []
WS_ERA = []

for file in filelist
    filenamesplit = split(file,"_")
    years = split(filenamesplit[end-1],"-")
    period = parse(Int,years[2])-parse(Int,years[1])
    if filenamesplit[4] == "ERAI" && filenamesplit[5] == "DD" && filenamesplit[end] =="beta-3.csv" && period > 30
        dir = CSV.read("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"*file,types=Float32,DataFrame)

        TIMESTAMP               = append!(TIMESTAMP,dir[!,"TIMESTAMP"][i:j])
        TA_ERA              = append!(TA_ERA,dir[!, "TA_ERA"][i:j])
        TA_ERA_NIGHT                = append!(TA_ERA_NIGHT,dir[!, "TA_ERA_NIGHT"][i:j])
        TA_ERA_NIGHT_SD                 = append!(TA_ERA_NIGHT_SD,dir[!, "TA_ERA_NIGHT_SD"][i:j])
        TA_ERA_DAY              = append!(TA_ERA_DAY,dir[!, "TA_ERA_DAY"][i:j])
        TA_ERA_DAY_SD               = append!(TA_ERA_DAY_SD,dir[!, "TA_ERA_DAY_SD"][i:j])
        SW_IN_ERA               = append!(SW_IN_ERA,dir[!, "SW_IN_ERA"][i:j])
        LW_IN_ERA               = append!(LW_IN_ERA,dir[!, "LW_IN_ERA"][i:j])
        VPD_ERA                 = append!(VPD_ERA,dir[!, "VPD_ERA"][i:j])
        PA_ERA              = append!(PA_ERA,dir[!, "PA_ERA"][i:j])
        P_ERA               = append!(P_ERA,dir[!, "P_ERA"][i:j])
        WS_ERA              = append!(WS_ERA,dir[!, "WS_ERA"][i:j])

    end
end

li = ["TIMESTAMP","TA_ERA","TA_ERA_NIGHT","TA_ERA_NIGHT_SD","TA_ERA_DAY",
"TA_ERA_DAY_SD","SW_IN_ERA","LW_IN_ERA","VPD_ERA","PA_ERA","P_ERA","WS_ERA"]
for (i,variable) in enumerate([TIMESTAMP,TA_ERA,TA_ERA_NIGHT,TA_ERA_NIGHT_SD,TA_ERA_DAY,
    TA_ERA_DAY_SD,SW_IN_ERA,LW_IN_ERA,VPD_ERA,PA_ERA,P_ERA,WS_ERA])
    variable = reshape(Float32.(variable),(N,Int(length(variable)/N)))

    jldsave(put_in*"F32_$(li[i]).jld2",
    data = variable)

end

"""
 "F32_LW_IN_ERA.jld2"
 "F32_PA_ERA.jld2"
 "F32_P_ERA.jld2"
 "F32_SW_IN_ERA.jld2"
 "F32_TA_ERA.jld2"
 "F32_TA_ERA_DAY.jld2"
 "F32_TA_ERA_DAY_SD.jld2"
 "F32_TA_ERA_NIGHT.jld2"
 "F32_TA_ERA_NIGHT_SD.jld2"
 "F32_TIMESTAMP.jld2"
 "F32_VPD_ERA.jld2"
 "F32_WS_ERA.jld2"

load("/net/scratch/lschulz/FLUXERAI/F32_TA_ERA_NIGHT.jld2")["data"]
"""

#in the fullset we have soil-temperature, gpp, nee and so on
#29 sites with length of longer then 15 years , need to be cut at 2005-2020
names = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")
i=0
for name in names[5:end-2]
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    years = parse(Int,filenamesplit[end-1][end-3:end]) - parse(Int,filenamesplit[end-1][end-8:end-5])
    y = years >= 15
    if  f && d && y
        println(name," $years")
        i+=1
    end
end

names = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")
i=0
for name in names[5:end-2]
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    metaindex = findall(x->x==String(filenamesplit[2]),meta[!,"SITE_ID"])
    years = parse(Int,filenamesplit[end-1][end-3:end]) - parse(Int,filenamesplit[end-1][end-8:end-5])
    y = years >= 15
    if  f && d && y && !isempty(metaindex)
        println(name," $years")
        println("class ",meta[metaindex,"IGBP_class"])
        println("height ",meta[metaindex,"height"][1])
        println("elevation ",meta[metaindex,"elevation"][1])
        i+=1
    end
end

"""
the spot variables
"""

#bsv <10
#cro crops
#csh >60
#cvm mixed crop
#dbf forest
#dnf forest
#ebf forest
#enf forest 
#gra grass
#mf forest
#osh 10-60
#sav 10-30
#sno ice
#urb urban
#wat water
#wet wetland
#wsa woody savannah

# some index of vegetation cover

["SNO","WAT","URB","BSV","WET","GRA","CVM","CRO","SAV","OSH","WSA","CSH","DBF","DNF","EBF","ENF","MF"]


"""
the choosing of observations
"""

#using this we proceed to write the data, lets first look at the available variables
s=CSV.read("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/FLX_DE-Tha_FLUXNET2015_FULLSET_DD_1996-2020_beta-3.csv",DataFrame)
s. + TTAB
lists all the variables

#this gives us a list of all the files with years longer then 15 and also available meta information
datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"
nameslist = readdir(datadir)
meta = CSV.read("logs/KW_38/meta_sites_unique.csv",DataFrame)
i=0
L = String[]
for name in nameslist[5:end-2]
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    metaindex = findall(x->x==String(filenamesplit[2]),meta[!,"SITE_ID"])
    years = parse(Int,filenamesplit[end-1][end-3:end]) - parse(Int,filenamesplit[end-1][end-8:end-5])
    y = years >= 15
    if  f && d && y && !isempty(metaindex)
        L = push!(L,name)
        i+=1
    end
end

#now we want to iterate over all these files and count the number of missing values in the specific column



for i in 1:size(s)[2]
    missings = String[]
    for spot in L
        series = CSV.read(datadir*spot,DataFrame)[:,i]
        if Int64(length(series)) - length(findall(x->x==-9999,series)) <= 5478
            missings = push!(missings,spot)
        end
    end
    println(i," ",missings)
end

#where do we have no missing values?
15:23 TA_F     TA_F_QC  TA_F_NIGHT  TA_F_NIGHT_SD  TA_F_NIGHT_QC  TA_F_DAY  TA_F_DAY_SD  TA_F_DAY_QC  SW_IN_POT
27:28 SW_IN_F  SW_IN_F_QC
32:33 LW_IN_F  LW_IN_F_QC
37:38 LW_IN_JSB_F  LW_IN_JSB_F_QC
42:43 VPD_F    VPD_F_QC (second probably broken)



#this gives us a filter of the list which reduces to complete variable set
datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"
names = readdir(datadir)
meta = CSV.read("logs/KW_38/meta_sites_unique.csv",DataFrame)
i=0

for variable in variable_list
    spots=0
    for spot in L[2:end]
        if isempty(findall(x->x==-9999,CSV.read(datadir*spot,DataFrame)[!,variable]))
            spots +=1
        else
            println(spot)
        end
    end
    println(variable," ",spots)
end


variable_list = ["SWC_F_MDS_1","SWC_F_MDS_1_QC","SWC_F_MDS_2","TS_F_MDS_1","TS_F_MDS_2","TS_F_MDS_3"]

variable_list = ["GPP_DT_VUT_MEAN","NEE_VUT_MEAN","RECO_DT_VUT_MEAN","SWC_F_MDS_1","TS_F_MDS_1"]

for (l,spot) in enumerate(L)
    print(split(spot,"_")[2]," ")
    csv = CSV.read(datadir*spot,DataFrame)
    filenamesplit = split(spot,"_")
    years = parse(Int,filenamesplit[end-1][end-3:end]) - parse(Int,filenamesplit[end-1][end-8:end-5])
    for variable in variable_list
        if !isempty(findall(x->x==variable,names(csv)))
            series = csv[!,variable]
            if length(series) - length(findall(x->x==-9999,series)) > 5478
                print(variable," ")
            end
        end
    end
    print(l," ",years,"\n")
end

4,5,18
L[deleteat!(Array(1:25),[4,5,18])]

"""
final program that gives us 15 years of data for our variables
one big matrix for the data
one small matrix for the meta data
    with the loop we seek th eindices of the meta file and gather the series
"""

# all 5 variables are not present in [1,4,5,9,11,18,19,25]
# oil temperature should at least be avilable in all but 4 5 13 18 24
indices = [4,5,13,18,24]
L_i = L[deleteat!(Array(1:25),indices)]
variable_list = ["GPP_DT_VUT_MEAN","NEE_VUT_MEAN","RECO_DT_VUT_MEAN","TS_F_MDS_1","TA_F","VPD_F"] #"SWC_F_MDS_1"
datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"
savedir = "/net/scratch/lschulz/fluxnetfullset/"

datafile = Array{Float32}(undef,5478,length(L_i),length(variable_list))
metaindices = Int64[]
for (l,name) in enumerate(L_i)
    csv = CSV.read(datadir*name,DataFrame)
    for (i,v) in enumerate(variable_list)
        datafile[:,l,i] = csv[end-5478+1:end,v]
    end
    filenamesplit = split(name,"_")
    f = filenamesplit[4] == "FULLSET"
    d = filenamesplit[5] == "DD"
    metaindices = push!(metaindices,findall(x->x==String(filenamesplit[2]),meta[!,"SITE_ID"])[1])
    #println(meta[metaindex,:])
    
end


findall(x->x==-9999,datafile[:,4,:]) #none  now!

sitesinfo = DataFrame(meta[metaindices,:])
wholedata = Array{Float32}(datafile)
jldsave(savedir*"fullset_15a_gpp_nee_reco_ts_ta_vpd.jld2",wholedata = wholedata,meta = sitesinfo)
# gives "/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2" wholedata.meta


"""
LAIF READOUT
"""

using Zarr

zz  = zopen("/net/data/LAI/LAI_AVHRR.zarr",fill_as_missing=true)
x = zz.arrays["layer"]

data = @view x[2101:2400,401:550,:]
datacols = reshape(data[:],size(data,1)*size(data,2),912)

function sectionlength(coli::SubArray{Union{Missing,Float32}})
    L = 0
    l = 0
    for i in coli
        l +=1
        L = (l>L) ? l : L
        if ismissing(i)
            l=0
        end
    end
    return L
end

l = [sectionlength(@view datacols[i,13:end]) for i=1:size(datacols,1)]

#with this we have 900 x 2 weeks of observation at the same time at 10.000 different spots...
sery = Float32.(datacols[findall(x-> x== 900,l),13:end])


l=  mapslices(x, dims = [1,2]) do i
    count(!ismissing, i)
end
f,ax,s = series(l[:]')
save(dir*"test.png",f)

"""
bbox-east-long 	180.0
bbox-north-lat 	90.0
bbox-south-lat 	-90.0
bbox-west-long 	-180.0
"""

#size 4320 longitude (degrees east) x 2160 latitude (degrees north) x 912

lon = Array(range(-180,180,4320))
lat = Array(range(90,-90,2160))

ger_latlot = [10.679,51.822]
ger_ind = [findall(x->abs(x-ger_latlot[1])<0.1,lon),findall(x->abs(x-ger_latlot[2])<0.1,lat)][1]

x[ger_ind[1],ger_ind[2],:]

#zarray takes much too long, lets use NetCDF

#one year
L = []
for year in 1981:2018
    x=ncread("/net/data/LAI/LAI_AVHRR/LAI.AVHRR.4320.2160.$(year).nc","LAI")
    #find the lat_long that have at least something
    nonempty = findall(x->x!=-9999.0,x)
    uno = unique([[nonempty[i][1],nonempty[i][2]] for i=1:length(nonempty)])
    L = append!(L,uno)
end

#convert to dataframe when there is this -9999 problem
df = DataFrame(Matrix(reshape(Array(x),size(x)[1]*size(x)[2],size(x)[3])'),:auto)
r = hcat([replace(x -> x==-9999.0 ? missing : x,df[:,k]) for k=1:size(df)[2]]...)