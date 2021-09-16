#' Characteristics of study participants.
#'
#' ART - antiretroviral therapy, CPT = co-trimoxazole preventative therapy
#'
#' @format A data frame with 425 rows and 32 variables:
#' \describe{
#'   \item{pid}{Unique patient identifier}
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
#'   \item{pid}{Unique patient identifier}
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


#' Results of stool testing for ESBL
#'
#' Results of testing stool/rectal swab samples for presence of any ESBL
#' producer
#'
#' @format A data frame with 1416 rows and 9 variables:
#' \describe{
#'   \item{pid}{Unique patient identifier}
#'   \item{lab_id}{Unique speciment identfier}
#'   \item{arm}{Study arm; 1 = sepsis, 2 = inpatient, no abx, 3 = community }
#'   \item{visit}{Study visit: 0 = baseline, 1 = d7, 2 = d28, 3 = d90, 4 = d180}
#'   \item{data_date}{Date of sample collection}
#'   \item{sample_type}{Sample type (stool or rectal swab}
#'   \item{ESBL}{Phenotypic ESBL presence (Positive) or absence (Negative)}
#'   \item{enroll_date}{Date of patient enrolment to study}
#'   \item{t}{Number of days post-enrolment that sample was collected}
#' }
"btESBL_stoolESBL"


#' Results of stool testing for ESBL
#'
#' Organisms identified from culture of stool on CHROMagar
#' (selective agar for cephalosporin resistance)
#'
#' @format A data frame with 1032 rows and 3 variables:
#' \describe{
#'   \item{lab_id}{Unique speciment identfier}
#'   \item{organism}{Organism identification}
#'   \item{ESBL}{Phenotypic presence of absence of ESBL in organism}
#' }
"btESBL_stoolorgs"

#' Posterior of fitted carriage model (decaying antimicrobial effect)
#'
#' Posterior of fitted longitudinal model of ESBL colonisation incorporating
#' an effect of antimicrobial exposure that exponentially decays following
#' cessation of exposure
#'
#' See manuscript for interpretation of parameters etc.
#'
#'
#' @format An S4 stanfit object (package "rstan") including model log-likelihood
#'
#'
"btESBL_model2posterior"

#' Posterior of fitted carriage model (constant antimicrobial effect)
#'
#' Posterior of fitted longitudinal model of ESBL colonisation incorporating
#' an effect of antimicrobial exposure that is stepwise-constant and
#' finsihes immediately following exposure
#'
#' See manuscript for interpretation of parameters etc.
#'
#'
#' @format An S4 stanfit object (package "rstan") including model log-likelihood
#'
#'
"btESBL_model1posterior"


#' Results of simulations from the model posterior
#'
#' These simulations set the antimicrobial and hospital exposure to
#' arbitrary values, set starting probability of colonisation to 0.5 then
#' draw parameter values from the posterior of the fitted model
#' (btESBL_model2posterior) to generate probability of being ESBL colonised
#' a time t later
#'
#'
#' @format A data frame with 600,000 rows and 18 variables:
#' \describe{
#'   \item{time}{Number of days post enrollment}
#'   \item{1}{Probability of being ESBL non-colonised at current time}
#'   \item{0}{Probability of being ESBL colonised at current time}
#'   \item{draw}{Unique ID of draw from posterior (1,2...1000)}
#'   \item{hosp_days}{Number of days of hospitalisation}
#'   \item{abx_days}{Number of days of antimicrobial exposure}
#'   \item{p0}{Probability at time 0 of being uncolonised}
#'   \item{p1}{Probability at time 0 of being colonised}
#'   \item{abx_start}{Time of antimicrobial exposure start}
#'   \item{abx_stop}{Time of antimicrobial exposure stop}
#'   \item{prev_abx}{Time of cessation of prev antimicrobial exposure (999=none)}
#'   \item{hosp_start}{Time of hospitalisation start}
#'   \item{hosp_stop}{Time of hospitalisation stop}
#'   \item{prev_hosp}{Time of cessation of prev hospitalisation (999=none)}
#' }
"btESBL_model2simulations"


#' Plasmid replicons identified in sequenced isolates
#'
#' Identified plasmid replicons using ARIBA and the PlamidFinder database.
#' See manuscript for details
#'
#' @format A tidy data frame with 2893 rows and 3 variables:
#' \describe{
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{ref_seq}{Identified plasmid replicon sequence}
#'   \item{species}{Species of sample}
#' }
"btESBL_plasmidreplicons"

#' Contig cluster membership
#'
#' Contig cluster assignment of all ESBL gene containing contigs, as
#' determined by cd-hit
#'
#' @format A tidy data frame with 714 rows and 10 variables:
#' \describe{
#'   \item{id}{Contig identifier}
#'   \item{clstr_size}{Number of contigs in cluster}
#'   \item{length}{Contig length, bases}
#'   \item{clstr_rep}{Is this the cd-hit cluster representative sample (1=yes,0=no)}
#'   \item{clstr_iden}{Identity (%) of sample to cluster representative sample}
#'   \item{clstr_cov}{Coverage (%) of sample to cluster representative sample}
#'   \item{gene}{ESBL gene}
#'   \item{clstr_name}{cd-hit cluster identfier: gene.n where n = 1,2..N and N is total no. clusters for given gene}
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{species}{Species of sample}
#' }
"btESBL_contigclusters"


#' Metadata of sequenced isolates
#'
#' Metadata of all sequenced isolates
#'
#' @format A tidy data frame with 676 rows and 12 variables:
#' \describe{
#'   \item{accession}{Sample NCBI accession number}
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{supplier_name}{Unique sample id}
#'   \item{pid}{Unique participant identifier}
#'   \item{arm}{Arm of study 1=sepsis, 2=inpatient, no abx, 3=community}
#'   \item{visit}{Study visit: 0 = baseline, 1 = d7, 2 = d28, 3 = d90, 4 = d180}
#'   \item{data_date}{Date of sample collection}
#'   \item{enroll_date}{Date of patient enrolment to study}
#'   \item{assess_type}{Number of days post-enrolment that sample was collected}
#'   \item{hosp_assoc}{Is sample community. hospital associated, or recent dc (see manuscript for details)}
#'   \item{hospoutcomedate}{Date of discharge from hospital. NA = never admittted}
#'   \item{Cluster}{PopPUNK cluster assignment. Prefix K = K. pneumo, E = E. coli}
#' }
"btESBL_sequence_sample_metadata"

#' E. coli core gene tree
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for E. coli
#'
#'
#' @format A "phylo" object ("ape" package)
#'
#'
"btESBL_coregene_tree_esco"


#' K. pneumoniae complex core gene tree
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for K. pneumoniae
#' complex
#'
#'
#' @format A "phylo" object ("ape" package)
#'
#'
"btESBL_coregene_tree_kleb"

#' Pairwise SNP-distance matrix for all E. coli genomes
#'
#' First column is sample names - all other columns are sample names
#'
#' @format A 473 x 474 tibble.
#'
#'
"btESBL_snpdists_esco"

#' Pairwise SNP-distance matrix for all K. pneumoniae complex genomes
#'
#' First column is sample names - all other columns are sample names
#'
#' @format A 350 x 351 tibble.
#'
"btESBL_snpdists_kleb"



