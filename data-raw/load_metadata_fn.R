# load metadata for genomes analysis




load_DASSIM3_metadata <- function(lanes_to_include, 
                          sanger_metadata_file,
                          location_of_phd_loading_script,
                          location_of_lims_loading_script,
                          return_which_df = "sanger",
                          recent_dc_def = 28) {
  
  # lanes to include as txt, one lane per line
  # sanger metadata file as .csv
  # location of scripts as string
  
#  require(tidyverse)
  
  source(location_of_phd_loading_script)
  source(location_of_lims_loading_script)
  
  # add in hospial associated vars
  
  lims_dates %>% 
    left_join(
      select(oc,pid,hospoutcomedate)
    ) ->
    lims_dates
  
  lims_dates$enroll_hospital_admission <- if_else(lims_dates$arm !=3,
                                                  lims_dates$enroll_date,
                                                  NA_Date_)
  
  lims_dates$enroll_hospital_discharge <- lims_dates$hospoutcomedate
  
  hosp_adms_fu <- 
    followup %>% 
    filter(d2visitnhospadmwhenadm != "") %>%  
    mutate(d2visitnhospaadmwhendisc =  
             dmy(gsub(" 00:00:00", "",d2visitnhospaadmwhendisc)),
           d2visitnhospadmwhenadm =  
             dmy(gsub(" 00:00:00",  "", d2visitnhospadmwhenadm)),
           d2visitnhospaadmwhendisc = 
             if_else(is.na(d2visitnhospaadmwhendisc),
                     data_date, 
                     d2visitnhospaadmwhendisc)) %>% 
    rename(fu_hosp_adm_date = d2visitnhospadmwhenadm,
           fu_dc_date = d2visitnhospaadmwhendisc) %>% 
    dplyr::select(pid, fu_hosp_adm_date, fu_dc_date)
  
  hosp_adms_enr <-
    enroll %>% 
    filter(pmhxrechospital == "Yes") %>% 
    dplyr::select(pid,
                  pmhxrechospital) %>% 
    rename(
      admitted_4wk_prior_enroll = pmhxrechospital)
  
  lims_dates %>% 
    left_join(
      hosp_adms_enr
    ) -> lims_dates
  
  lims_dates %>% 
    mutate(admitted_4wk_prior_enroll = 
             if_else(is.na(admitted_4wk_prior_enroll),
                     "No",
                     admitted_4wk_prior_enroll)) ->
    lims_dates
  
  
  lims_dates %>% 
    left_join(
      hosp_adms_fu
    ) -> lims_dates
  
  # if there are two dates matched to one lab id - pick the closest that is before
  # the data date
  # ie hosp_adm needs to be before data_date
  # and then the closest one
  # so make two vars
  # one - code o for adm date before data_date, 1 otherwise
  # then another abs(diff)
  # arrnge by (before_var, abd(dff)) %>%  slice(n=1)
  
  lims_dates %>% mutate(
    hosp_adm_after_data_date = as.numeric(fu_hosp_adm_date > data_date),
    abs_diff_hosp_adm_date_data_date = abs(fu_hosp_adm_date - data_date)) %>% 
    group_by(lab_id) %>% 
    arrange(lab_id,hosp_adm_after_data_date, 
            abs_diff_hosp_adm_date_data_date) %>% 
    slice(n=1) -> lims_dates
  
  
  # now - 
  # isolate is IN HOSPITAL if 
  # data_date > enroll date AND data_date <= discharge_date AND arm %in% c(1,2)
  # OR
  # data_date > fu_hosp_adm_date AND data_date <= fu_dc_date
  
  # isolate is - RECENT DC if
  # data_date - enroll discharge date < 4 week and 
  # OR
  # data_date - fu dischareg date < 4 week
  # OR
  # admitted_4wk_prior_enroll == "Yes and enrollment sample
  
  # otherwise community
  
  lims_dates %>% 
    mutate(hosp_assoc =
             case_when(
               data_date > enroll_hospital_admission &
                 data_date <= enroll_hospital_discharge &
                 arm %in% c(1,2,4) ~ "in_hospital",
               data_date > fu_hosp_adm_date &
                 data_date <= fu_dc_date ~ "in_hospital",
               
               visit == 0 & 
                 admitted_4wk_prior_enroll == "Yes" ~ "recent_dc",
               as.numeric(data_date - 
                            enroll_hospital_discharge) <= recent_dc_def &
                 as.numeric(data_date - 
                              enroll_hospital_discharge) > 0 ~ "recent_dc",
               as.numeric(data_date - 
                            fu_dc_date) <= recent_dc_def &
                 as.numeric(data_date - 
                              fu_dc_date) > 0 ~ "recent_dc",
               
               TRUE ~ "community")
    ) -> lims_dates
  
  
  
  lanes <- read_lines(lanes_to_include)
  sang.metadata <- 
    read_csv(sanger_metadata_file) %>% 
    rename_with(~ tolower(gsub(" ", "_", .x))) %>% 
    mutate(supplier_name = str_split(
      supplier_name, "_", simplify = TRUE)[,1]) 
  
 
  
  # Make a unified follow up df by combinging enroll and followup visits
  
  sang.metadata <-
    filter(sang.metadata, lane %in% lanes)

  left_join(sang.metadata, select(lims_dates, 
                                  pid, arm, lab_id, 
                                  visit, 
                                  data_date, 
                                  enroll_date, 
                                  assess_type,
                                  hosp_assoc),
            by = c("supplier_name" = "lab_id")) -> 
    dfout 
if (return_which_df =="sanger") {
  return(dfout)
} else {
  return(lims_dates)
}
 }
