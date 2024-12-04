# extract relevant DHS/MIS output from all DHS downloaded files



rm(list = ls()) ##clean your environment (important so that you don't accidentally bind unrelated dataframes later)

require(readstata13)
require(plyr)
require(dplyr)
require(gdata)

setwd("C:/Users/moniqueam/Downloads/DHS_20220926")

identify_vars = c("v000", "hv000", #country-survey code
                  "v001", "hv001", "DHSCLUST", #cluster
                  "v005", #ind weights
                  "hv007", "v007",  # year
                  'v006', 'hv006',    # month of survey in cluster
                  "v012", "hv105" #respondent age
)
malaria_vars = c(
             # PfPR (microscopy)
             'hml32', 'hml32a', 'hml33', 'hml35', 
             # ITNs
             'hml20','hml19','s508', '460','461',
             # treatment-seeking
             'h32z',
             # receive heel prick given seek treatment
             'h47',
             # IPTp >=1 dose
             'm49a',
             # IPTp number of doses
             'ml1',
             # time between fever onset and ACT
             'ml20a',
             # child took combination with artemisinin
             #BASE for ML13A-Z: Children suffering from fever or short rapid breaths or difficulty
             # breathing in the last two weeks (H22 = 1 or H31B = 1). In previous recodes the base was
             # restricted to children with fever in the last 2 weeks (H22 = 1).
             # Questions pertaining to ML14A to ML14Z are no longer part of the DHS VII core
             # questionnaire, but the variables are kept in the DHS VII recode
             'ml13e',
             # months ago household obtained net
             'hml4'
)

varlist = c(identify_vars, malaria_vars)


temp = list.files(pattern="*.DTA") #list all files in the working directory with extension

for (i in 1:length(temp)) assign(temp[i], read.dta13(temp[i], select.cols = varlist))  #bring in all files in folder with selected variables

dfs = sapply(.GlobalEnv, is.data.frame) #list all dataframes in the environment

dhs <- do.call(rbind.fill, mget(names(dfs)[dfs])) #bind all the dataframes in environment, missing columns filled in

keep(dhs, sure = T) #remove all datasets except dhs, to preserve memory/space

# dhs$country_code <- gsub("[[:digit:]]","",dhs$v000)
# dhs$survey_round <- gsub("[^[:digit:]]", "", dhs$v000)

# divide into household versus individual
dhs_ind = dhs[!is.na(dhs$v001),]
dhs_hh = dhs[!is.na(dhs$hv001),]

# remove rows that don't have any of the malaria information needed (i.e., all NA for variables related to malaria or interventions)
# individual df
malaria_vars_df = intersect(malaria_vars, colnames(dhs_ind))
dhs_ind2 = dhs_ind[rowSums(is.na(dhs_ind[,which(colnames(dhs_ind) %in% malaria_vars)])) != length(malaria_vars_df), ]
dim(dhs_ind2)
View(dhs_ind2)
# household df
malaria_vars_df = intersect(malaria_vars, colnames(dhs_hh))
dhs_hh2 = dhs_hh[rowSums(is.na(dhs_hh[,which(colnames(dhs_hh) %in% malaria_vars)])) != length(malaria_vars_df), ]
dim(dhs_hh2)
View(dhs_hh2)



# subset to sites near study locations




save(dhs, file = "dhs.csv")














##################################################################################
# determine which clusters are in which chiefdoms
##################################################################################


# turn MIS_2016 into spatial points data frame
# not sure what the spatial projection is... will try with same spatial projection as shapefile
points_crs = crs(admin_shape)
MIS_2016_shape = SpatialPointsDataFrame(MIS_2016[,c('longitude', 'latitude')],
                                        MIS_2016,
                                        proj4string = points_crs)
# find which chiefdom each cluster belongs to
MIS_2016_locations = over(MIS_2016_shape, admin_shape)
if(nrow(MIS_2016_locations) == nrow(MIS_2016_shape)){
  MIS_2016_shape$NAME_3 = MIS_2016_locations$NAME_3
  MIS_2016_shape$NAME_2 = MIS_2016_locations$NAME_2
  MIS_2016_shape$NAME_1 = MIS_2016_locations$NAME_1
}

write.csv(as.data.frame(MIS_2016_shape), 'C:/Users/mambrose/Dropbox (IDM)/Malaria Team Folder/projects/SierraLeone_hbhi/explore_DHS/MIS_2016.csv')


par(mfrow=c(3,4))
# for each of the interventions/measures, count number of individuals surveyed in each chiefdom
num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_3','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
  group_by(NAME_3) %>%
  summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
            num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
            num_w_fever=sum(num_w_fever, na.rm = TRUE),
            num_preg=sum(num_preg, na.rm = TRUE))
num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))


# look at survey numbers by district
# for each of the interventions/measures, count number of individuals surveyed in each chiefdom
num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_2','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
  group_by(NAME_2) %>%
  summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
            num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
            num_w_fever=sum(num_w_fever, na.rm = TRUE),
            num_preg=sum(num_preg, na.rm = TRUE))
num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,800, length.out=80))
hist(num_surveyed$num_use_itn_all,  main='ITN', breaks=seq(0,1600, length.out=80))
hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,200, length.out=80))
hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,500, length.out=80))

# look at survey numbers by region
# for each of the interventions/measures, count number of individuals surveyed in each chiefdom
num_surveyed = as.data.frame(MIS_2016_shape[c('NAME_1','num_mic_test', 'num_use_itn_all', 'num_w_fever', 'num_preg')]) %>% 
  group_by(NAME_1) %>%
  summarise(num_mic_test=sum(num_mic_test, na.rm = TRUE),
            num_use_itn_all=sum(num_use_itn_all, na.rm = TRUE),
            num_w_fever=sum(num_w_fever, na.rm = TRUE),
            num_preg=sum(num_preg, na.rm = TRUE))
num_surveyed = num_surveyed[!is.na(num_surveyed[,1]),]
hist(num_surveyed$num_mic_test,main='microscopy (PfPR)', breaks=seq(0,3000, length.out=80))
hist(num_surveyed$num_use_itn_all,  main='ITN',breaks=seq(0,6000, length.out=80))
hist(num_surveyed$num_w_fever, main='fever (CM)', breaks=seq(0,700, length.out=80))
hist(num_surveyed$num_preg, main='preganacies (IPTp)', breaks=seq(0,2000, length.out=80))
par(mfrow=c(1,1))




# get weighted means and number tested/positive within each chiefdom for all variables



