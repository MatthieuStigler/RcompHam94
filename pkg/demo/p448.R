library(fracdiff)
args(fdGPH)


data(gnptbill,package = "RcompHam94")
print( fdGPH(gnptbill$GNP) )
print( fdGPH(gnptbill$TBILL) )


