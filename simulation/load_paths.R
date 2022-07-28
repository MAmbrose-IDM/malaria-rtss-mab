get_project_paths = function(){
  if (!exists("user_path")) user_path <- gsub("Documents", "", Sys.getenv("HOME"))
  cond1 <- length(grep("b1139", user_path)) > 0
  cond2 <- length(grep("home", user_path)) > 0
  
  if (cond1 | cond2) {
    projectpath <- '/projects/b1139/rtss-scenario_IO/'
    datapath <- file.path(projectpath, 'data')
  }else {
    home_path <- file.path(user_path, 'Box', 'NU-malaria-team')
    datapath <- file.path(home_path, 'data')
    projectpath <- file.path(home_path, 'projects', 'rtss_scenarios')
    if (!dir.exists(projectpath)) {
      projectpath <- file.path(user_path, 'Box', 'rtss_scenarios')
    }
    if(grepl('moniqueam', user_path)){
      home_path = "C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/projects/mAb_rtss_smc_comparison"
      datapath = file.path(home_path, 'reference_data')
      projectpath = home_path
    }
  }
  return(c(datapath, projectpath))
}