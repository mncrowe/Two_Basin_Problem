using NetCDF

filename = "QG.nc"

q = ncread(filename, "Q")
psi = ncread(filename, "psi")
E = ncread(filename, "E")
Z = ncread(filename, "Z")

nothing
