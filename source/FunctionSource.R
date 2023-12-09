 
setwd(paste(directory,"/source/", sep = ''))    # set temp working directory 

source(paste(getwd(), "/Replicates.R", sep = ''))
source(paste(getwd(), "/BetaParams.R", sep = ''))
source(paste(getwd(), "/max1.R", sep = ''))

setwd(directory) #change directory back