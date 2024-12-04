# only keep simulation sweep4c that weren't already in 4 or 4b

require(sqldf)

# project filepath
projectpath = 'C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/projects/mAb_rtss_smc_comparison'

coordinator_4 = read.csv(paste0(projectpath, '/simulation_inputs/coordinator_files/sweep4.csv'))
coordinator_4b = read.csv(paste0(projectpath, '/simulation_inputs/coordinator_files/sweep4b.csv'))
coordinator_4c = read.csv(paste0(projectpath, '/simulation_inputs/coordinator_files/sweep4c.csv'))

coordinator_4 = coordinator_4[,-which(colnames(coordinator_4) =='setting_id')]
coordinator_4b = coordinator_4b[,-which(colnames(coordinator_4b) =='setting_id')]
coordinator_4c = coordinator_4c[,-which(colnames(coordinator_4c) =='setting_id')]

# total scenarios already run
dim(coordinator_4)[1] + dim(coordinator_4b)[1]
# total scenarios in the full eventual run (4c)
dim(coordinator_4c)[1]

# need to remove the scenarios already run from 4c to get the new scenarios that need to run
coordinator_4d = sqldf('SELECT * FROM coordinator_4c EXCEPT SELECT * FROM coordinator_4')
coordinator_4d = sqldf('SELECT * FROM coordinator_4d EXCEPT SELECT * FROM coordinator_4b')
dim(coordinator_4d)[1]


# add the setting_id column
coordinator_4d$setting_id = paste0('HX', seq(0,nrow(coordinator_4d)-1))
write.csv(coordinator_4d, paste0(projectpath, '/simulation_inputs/coordinator_files/sweep4d.csv'), row.names = FALSE)
