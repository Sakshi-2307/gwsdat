library(GWSDAT)
opt <- createOptions("Site Name")
opt$WellDataFilename <- '/Users/sakshiprasad/Downloads/GWSDAT Online Input Data sets/BasicExample_WellData.csv'
opt$WellCoordsFilename <- '/Users/sakshiprasad/Downloads/GWSDAT Online Input Data sets/BasicExample_WellCoords.csv'
launchApp(opt)
