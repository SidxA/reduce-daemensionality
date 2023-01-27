#time stamp for saving analysis results
using Dates


function init_logging()
    weekno = week(unix2datetime(time()))
    datestring = string("KW_",lpad(weekno,2,"0"),"/")
    workdir = "/net/home/lschulz/logs/"
    dir = workdir*datestring
    if isdir(dir)==false
        mkdir(dir)
    end
    return dir
end