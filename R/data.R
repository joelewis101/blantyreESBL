R#' Characteristics of study participants.
#'
#' ART - antiretroviral therapy, CPT = co-trimoxazole preventative therapy
#'
#' @format A data frame with 425 rows and 32 variables:
#' \describe{
#'   \item{pid}{Unique patient identified}
#'   \item{arm}{Study arm; 1=sepsis, 2=inpatient, no abx, 3=community}
#'   \item{enroll_date}{Date of enrolment in study}
#'   \item{calc_age}{Age, years, at enrolment}
#'   \item{ptsex}{Self reported gender}
#'   \item{hivstatus}{HIV status}
#'   \item{recieved_prehosp_ab}{Received antimicrobials (except CPT) last 4wk?}
#'   \item{pmhxrechospital}{Hospital admission last 4wk?}
#'   \item{tbstatus}{Ever been treated for TB?}
#'   \item{tbongoing}{If tbstatus=Yes, currently being treated for TB?}
#'   \item{hivonart}{If hivstatus=Reactive, current taking ART?}
#'   \item{hivart}{If currently taking ART, which regimen?}
#'   \item{art_time}{If currently taking ART, months on ART}
#'   \item{hivcpt}{If hivstatus=Reactive, currently taking CPT?}
#'   \item{tobaccoyn}{Currently smokes tobacco}
#'   \item{alcoholyn}{Currently drinks alcohol}
#'   \item{highestedu}{Highest education attained}
#'   \item{job}{Current employment, if any}
#'   \item{householdadultsno}{Number of adults in participant's household}
#'   \item{householdchildno}{Number of children in participant's household}
#'   \item{toilet}{Type of toilet in participant's household}
#'   \item{watersource}{Usual household water source for drinking}
#'   \item{watertreated}{Is drinking water in the household usually treated?}
#'   \item{electricityyn}{Electricity available in household?}
#'   \item{fuel}{Usual cooking fuel}
#'   \item{keepanim}{Are animals kept in the household?}
#'   \item{keep.poultry}{If animals are kept, are poultry kept?}
#'   \item{keep.goats}{If animals are kept, are goats kept?}
#'   \item{keep.dogs}{If animals are kept, are dogs kept?}
#'   \item{keep.cattle}{If animals are kept, are cattle kept?}
#'   \item{keep.sheep}{If animals are kept, are sheep kept?}
#'   \item{keep.mules}{If animals are kept, are mules kept?}
#' }
"btESBL_participants"


#' Exposures of included participants
#'
#' These data are recorded in long format. Covariate values are included
#' in columns, and there is one row for each day where a covariate value
#' changes or when the participant has their final observation (i.e date
#' of censoring). So, if a participant has a covariate value that changes
#' from 1 to 0 on subsequent rows, the covariate value should be
#' interpreted as being 1 up until the day of change
#'
#' @format A data frame with 1838 rows and 26 variables:
#' \describe{
#'   \item{pid}{Unique patient identified}
#'   \item{assess_type}{Day of observation (day 0=enrolment day)}
#'   \item{died}{1=died on day of current observation}
#'   \item{hosp}{1=hospitalised on day of current observation, 0 =not}
#'   \item{amoxy}{Exposure to amoxicillin on day of observation}
#'   \item{genta}{Exposure to gentamicin on day of observation}
#'   \item{azithro}{Exposure to azithromycin on day of observation}
#'   \item{tb}{Exposure to first-line TB treatment on day of observation}
#'   \item{benzy}{Exposure to penicillin on day of observation}
#'   \item{fluco}{Exposure to fluconazole on day of observation}
#'   \item{ampho}{Exposure to amphotericin on day of observation}
#'   \item{quin}{Exposure to quinine on day of observation}
#'   \item{chlora}{Exposure to chloraquine on day of observation}
#'   \item{coart}{Exposure to artemether/lumefantrine on day of observation}
#'   \item{cipro}{Exposure to ciprofloxacin on day of observation}
#'   \item{coamo}{Exposure to co-amoxiclav on day of observation}
#'   \item{arte}{Exposure to artesunate on day of observation}
#'   \item{cotri}{Exposure to co-trimoxazole on day of observation}
#'   \item{clinda}{Exposure to clindamicin on day of observation}
#'   \item{doxy}{Exposure to doxycycline on day of observation}
#'   \item{erythro}{Exposure to erythomycin on day of observation}
#'   \item{acicl}{Exposure to aciclovir on day of observation}
#'   \item{fluclox}{Exposure to flucloxacillin on day of observation}
#'   \item{metro}{Exposure to metronidazole on day of observation}
#'   \item{strepto}{Exposure to streptomicin on day of observation}
#' }
"btESBL_exposures"
