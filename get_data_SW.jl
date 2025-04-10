using NCDatasets

filename = "shallow_water.nc"

SWData = NCDataset(filename)

ω = (SWData["ω"][1:1024, 1:1024, :] .+ SWData["ω"][2:1025, 1:1024, :] .+
     SWData["ω"][1:1024, 2:1025, :] .+ SWData["ω"][1:1024, 2:1025, :]) / 4

close(SWData)