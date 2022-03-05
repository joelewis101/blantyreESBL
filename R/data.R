#' Data: Characteristics of study participants.
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


#' Data: Exposures of included participants
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


#' Data: Results of stool culture for ESBL
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


#' Data: Species of bacteria isolated from stool culture for ESBL
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

#' Data: Posterior of fitted carriage model (decaying antimicrobial effect)
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

#' Data: Posterior of fitted carriage model (constant antimicrobial effect)
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


#' Data: Posterior of fitted carriage model (seperate ceftriaxone effect)
#'
#' Posterior of fitted longitudinal model of ESBL colonisation incorporating
#' an effect of antimicrobial exposure that exponentially decays following
#' cessation of exposure, with ceftriaxone effect seperate from
#' all otehr antimicrobials
#'
#' See manuscript for interpretation of parameters etc.
#'
#'
#' @format An S4 stanfit object (package "rstan") including model log-likelihood
#'
#'
"btESBL_model3posterior"

#' Data: Results of simulations from the model posterior
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
#'   \item{sim_run}{Simulation run ID}
#'   \item{draw}{Unique ID of draw from posterior (1,2...1000)}
#'   \item{time}{Time t in days}
#'   \item{hosp_days}{Number of days of hospitalisation}
#'   \item{abx_days}{Number of days of antimicrobial exposure}
#'   \item{pr_esbl_pos_t0esblneg}{Probability of ESBL colonisation at time t conditional on no colonisation at time 0}
#'   \item{pr_esbl_pos_t0esblpos}{Probability of ESBL colonisation at time t conditional on colonisation at time 0}
#'   \item{pr_esbl_pos}{Overall estimated prevalence of ESBL colonisation at time t assuming initial prevalence 0.5}
#'   \item{prev_hosp}{Time of cessation of prev hospitalisation (999=none)}
#' }
"btESBL_model2simulations"


#' Data: Plasmid replicons identified in sequenced isolates
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

#' Data: AMR genes identified in sequenced isolates (excluding QRDR mutations)
#'
#' Identified AMR genes using ARIBA and the SRST database.
#' See manuscript for details
#'
#' @format A tidy data frame with 10822 rows and 3 variables:
#' \describe{
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{ref_seq}{Identified AMR gene sequence}
#'   \item{genus}{Genus of sample (E. coli or K. pneumioniae complex)}
#' }
"btESBL_amrgenes"

#' Data: Identified QRDR mutations
#'
#' Identified quinolone resistance determinant region
#'mutations in the samples of this study
#'
#' @format A tidy data frame with 23 rows and 3 variables:
#' \describe{
#'   \item{gene}{QDRR gene (gyrA, gyrB, parC, parE)}
#'   \item{variant}{Identified variant}
#'   \item{pmid}{Pubmed ID of publication describing mutation}
#' }
"btESBL_qrdr_mutations"

#' Data: CARD described QRDR mutations
#'
#' Quinolone resistnce determinant region mutations that have been
#' associated with resistance in E. coli, downloaded from the
#' Comprehensive Antibiotic Resistance Database https://card.mcmaster.ca/
#' on 6 September 2021
#'
#' @format A tidy data frame with 1293 rows and 4 variables:
#' \describe{
#'   \item{gene}{QDRR gene (gyrA, gyrB, parC, parE)}
#'   \item{variant}{Identified variant}
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{genus}{Genus of sample (E. coli or K. pneumioniae complex)}
#' }
"btESBL_CARD_qrdr_mutations"

#' Data: NCBI phenotypic beta-lactamase classifications
#'
#' Classification of beta lactamase genes downloaded from NCBI at
#' https://ftp.ncbi.nlm.nih.gov/pathogen/betalactamases/Allele.tab
#' on 17 September 2021
#'
#' @format A tidy data frame with 1293 rows and 4 variables:
#' \describe{
#'   \item{allele_name}{Allele name}
#'   \item{protein_accession_}{Protein accession}
#'   \item{nucleotide_accession_}{nucleotide accession}
#'   \item{gene_name}{Gene name}
#'   \item{curated_gene_product_name}{Curated phenotype}
#'   \item{class}{Inferred class of beta lactamase from phenotype}
#' }
"btESBL_NCBI_phenotypic_bl"

#' Data: Contig cluster membership
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


#' Data: Accession numbers and metadata of isolates sequenced for this study
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
#'   \item{number_of_contigs}{Number of contigs in assembly}
#'   \item{N50}{N50 of assembly}
#'   \item{ecoli_phylogroup}{E. coli phylogroup from in silico PCR}
#'   \item{ecoli_pathogroup}{E. coli pathogroup (see manuscript for details)}
#'   \item{species}{Isolate species either from API (E. coli) or Kleborate}
#'   \item{kleb_k_locus}{K. pnemo complex inferred K locus from Kleborate}
#'   \item{kleb_k_locus_confidence}{Kleborate defined K locus confidence}
#'   \item{kleb_o_locus}{K. pnemo complex inferred O locus from Kleborate}
#'   \item{kleb_o_locus_confidence}{Kleborate defined O locus confidence}
#'   \item{kleb_YbST}{Klebsiella YbST allele as per Kleborate}
#'   \item{kleb_CbST}{Klebsiella CbST allele as per Kleborate}
#'   \item{kleb_AbST}{Klebsiella AbST allele as per Kleborate}
#'   \item{kleb_SmST}{Klebsiella SmST allele as per Kleborate}
#'   \item{kleb_rmpA}{Klebsiella rmpA allele as per Kleborate}
#'   \item{kleb_rmpA2}{Klebsiella rmpA2 allele as per Kleborate}
#' }
"btESBL_sequence_sample_metadata"

#' Data: E. coli core gene tree
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for E. coli
#'
#'
#' @format A "phylo" object ("ape" package)
#'
#'
"btESBL_coregene_tree_esco"


#' Data: K. pneumoniae complex core gene tree
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for K. pneumoniae
#' complex
#'
#'
#' @format A "phylo" object ("ape" package)
#'
#'
"btESBL_coregene_tree_kleb"

#' Data: Pairwise SNP-distance matrix for all E. coli genomes
#'
#' First column is sample names - all other columns are sample names
#'
#' @format A 473 x 474 tibble.
#'
#'
"btESBL_snpdists_esco"

#' Data: Pairwise SNP-distance matrix for all K. pneumoniae complex genomes
#'
#' First column is sample names - all other columns are sample names
#'
#' @format A 350 x 351 tibble.
#'
"btESBL_snpdists_kleb"


#' Data: metadata of 10,146 global E. coli genomes from Horesh et al
#'
#' Metadata of a global collection of 10,146 E. coli genomes
#' Taken from Horesh et al 2021
#' See original Horesh publication for details and variable definition
#' https://doi.org/10.1099/mgen.0.000499
#'
#' @format A data frame with 10,146 rows and 22 variables:
#'
"btESBL_ecoli_horesh_metadata"

#' Data: Accession numbers and metadata of 97 E. coli genomes from Malawi from Musciha et al
#'
#' Metadata of a collection of 97 E. coli genomes selected for temporal
#' diversity and diversity in antibiogram
#' Taken from Musicha et al 2017
#' For full details see original Musicha publication
#' https://doi.org/10.1093/jac/dkx058
#'
#' The phylogroup and ST here were determined for this study and were not
#' provided in the original metadata; the popPUNK cluster assignemt uses the
#' Horesh database
#'
#' #' \describe{
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{ID}{Unmique sample ID}
#'   \item{Year}{Year of sample collection}
#'   \item{Sample_type}{Type of clincial sample from which isolate originates}
#'   \item{Age_group}{Age of patient from which sample orignated}
#'   \item{Organism}{Organism}
#'   \item{ST}{Multlocus sequence type as determined by ARIBA}
#'   \item{phylogroup}{Phylogroup following Clermon scheme with in silico PCR}
#'   \item{Cluster}{PopPUNK cluster assignment using the Horesh E. coli database}
#'   \item{Accession}{NCBI Accession number}
#'
#' }
"btESBL_ecoli_musicha_metadata"

#' Data: popPUNK cluster assignment of study E. coli isolates using global database
#'
#' PopPUNK cluster assignment of study E. coli using the popPUNK database
#' derived using a global E. coli collection as described in Horesh et al
#' 2021; see manuscript for details
#'
#' @format A tidy data frame with 473 rows and 2 variables:
#' \describe{
#'   \item{lane}{Unique sample-sequencing run ID}
#'   \item{Cluster}{popPUNK cluster assignment}
#' }
"btESBL_ecoli_global_popPUNK_clusters"

#' Data: E.coli global context tree
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for E. coli
#' incorporating study and global context isolates (see manuscript for
#' details)
#'
#' @format A "phylo" object ("ape" package)
#'
#'
"btESBL_ecoli_globaltree"

#' Data: K. pneumoniae complex global context phylogeny
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for K.
#' pneumoniae complex incorporating study and global context isolates
#'  (see manuscript for details)
#'
#' @format A "phylo" object ("ape" package)
#'
"btESBL_kleb_globaltree"


#' Data: Accession numbers and metadata of 682 global K. pneumoniae complex genomes
#'
#' Metadata of a global collection of 682 K. pneumoniae complex genomes.
#' 66 genomes from kenya from the Henson et al
#' https://doi.org/10.1016/j.ijmm.2017.07.006
#'
#' 72 genomes fro Malawi described in Musicha et al
#' https://doi.org/10.1093/jac/dkz032
#'
#' 100 genomies from Malawi as described in Cornick et al
#' https://doi.org/10.1101/2020.08.06.236117
#'
#' 288 global genomes from Holt et al
#'https://doi.org/10.1073/pnas.1501049112
#'
#'@format A dataframe with 682 rows and 23 variables
#' \describe{
#'   \item{name}{unique identifier}
#'   \item{ESBL}{ESBL gene present? ESBL= yes, 0 = no}
#'   \item{Sample accession}{NCBI accession number}
#'   \item{study}{study sample was taken from}
#'   \item{Sample Source}{Source of sample}
#'   \item{Isolate Type}{For human samples, are they carriage or infecting?}
#'   \item{location}{Lcoation (global, Kenyan, Malawian)}
#'   \item{ST}{Multilocus sequence type as determined by ARIBA}
#'   \item{K_locus}{K. pnemo complex inferred K locus from Kleborate}
#'   \item{K_locus_confidence}{Kleborate defined K locus confidence}
#'   \item{O_locus}{K. pnemo complex inferred O locus from Kleborate}
#'   \item{O_locus_confidence}{Kleborate defined O locus confidence}
#'   \item{YbST}{Klebsiella YbST allele as per Kleborate}
#'   \item{CbST}{Klebsiella CbST allele as per Kleborate}
#'   \item{AbST}{Klebsiella AbST allele as per Kleborate}
#'   \item{SmST}{Klebsiella SmST allele as per Kleborate}
#'   \item{rmpA}{Klebsiella rmpA allele as per Kleborate}
#'   \item{rmpA2}{Klebsiella rmpA2 allele as per Kleborate}
#'   \item{ybt}{Presence of yersiniabactin virulence locus (0=absent,1=present)}
#'   \item{clb}{Presence of colibactin virulence locus (0=absent,1=present)}
#'   \item{iuc}{Presence of aerobactin virulence locus (0=absent,1=present)}
#'   \item{iro}{Presence of salmochelin virulence locus (0=absent,1=present)}
#'   \item{Malawi}{Was sample collected in Malawi (0=no,1=yes)}
#'
#' }
"btESBL_kleb_global_metadata"

#' Data: K. pneumoniae complex core gene tree of Malawian carriage and invasive isolates
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for K.
#' pneumoniae complex incorporating all Malawian isolates
#'  (see manuscript for details)
#'
#' @format A "phylo" object ("ape" package)
"btESBL_kleb_malawi_allisolate_core_gene_tree"

#' Data: E. coli ST167 global phylogeny
#'
#' Midpoint rooted maximum likelihood core gene phylogeny for E. coli
#' incorporating all ST167 study isolates and all ST167 E. coli
#' downloaded from Enterobase on 1st March 2021 using multiple
#' sequence alignment from mapping to reference.
#'  (see manuscript for details)
#'
#' @format A "phylo" object ("ape" package)
"btESBL_ecoli_globalst167_tree"

#' Data: E. coli ST410 global phylogeny
#'
#' Midpoint rooted maximum likelihood phylogeny for E. coli
#' incorporating all ST410 study isolates and all ST410 E. coli
#' downloaded from Enterobase on 1st March 2021 using multiple
#' sequence alignment from mapping to reference.
#'  (see manuscript for details)
#'
#' @format A "phylo" object ("ape" package)
"btESBL_ecoli_globalst410_tree"

#' Data: E. coli ST131 global phylogeny
#'
#' Midpoint rooted maximum likelihood phylogeny for E. coli
#' incorporating all ST131 study isolates and all ST131 E. coli
#' from a study of global ST131
#' (Mcnally et al 2019 https://doi.org/10.1128/mBio.00644-19)
#' reconstructed from a multiple sequence alignment from mapping to reference.
#'  (see manuscript for details)
#'
#' @format A "phylo" object ("ape" package)
"btESBL_ecoli_globalst131_tree"

#' Data: metadata of global E. coli ST410 isolates
#'
#' Metadata of a collection of 1230 E. coli ST410 genomes
#' downloaded from Enterobase on 1st March 2021
#'@format A dataframe with 1230 rows and 12 variables
#'\describe{
#'   \item{Uberstrain}{Uberstrain as defined by Enterobase}
#'   \item{Name}{Sample ID}
#'   \item{accession}{Sample accession}
#'   \item{platform}{Sequencing platform}
#'   \item{library}{library as defined by Enterobase}
#'   \item{insert_size}{insert size}
#'   \item{Experiment}{Experiment accession number}
#'   \item{Source Niche}{Source niche as defined by Enterobase e.g. human , animal}
#'   \item{Source Details}{Source details as defined by Enterobase e.g. species}
#'   \item{Country}{Country of sample collection}
#'   \item{Collection Year}{Year of sample colllection}
#'   \item{ST}{Multilocus sequence type}
#' }
"btESBL_ecoli_st410_metadata"


#' Data: metadata of global E. coli ST167 isolates
#'
#' Metadata of a collection of 1230 E. coli ST167 genomes
#' downloaded from Enterobase on 1st March 2021
#'@format A dataframe with 1230 rows and 12 variables
#'\describe{
#'   \item{Uberstrain}{Uberstrain as defined by Enterobase}
#'   \item{Name}{Sample ID}
#'   \item{accession}{Sample accession}
#'   \item{platform}{Sequencing platform}
#'   \item{library}{library as defined by Enterobase}
#'   \item{insert_size}{insert size}
#'   \item{Experiment}{Experiment accession number}
#'   \item{Source Niche}{Source niche as defined by Enterobase e.g. human , animal}
#'   \item{Source Details}{Source details as defined by Enterobase e.g. species}
#'   \item{Country}{Country of sample collection}
#'   \item{Collection Year}{Year of sample colllection}
#'   \item{ST}{Multilocus sequence type}
#' }
"btESBL_ecoli_st167_metadata"


#' Data: AMR genes identified in global ST167 isolates
#'
#' Identified AMR genes using ARIBA and the SRST database.
#' See manuscript for details
#'
#' @format A tidy data frame with 5494 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{gene}{Identified AMR gene sequence}
#' }
"btESBL_ecoli_st167_amr"


#' Data: AMR genes identified in global ST410 isolates
#'
#' Identified AMR genes using ARIBA and the SRST database.
#' See manuscript for details
#'
#' @format A tidy data frame with 7704 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{gene}{Identified AMR gene sequence}
#' }
"btESBL_ecoli_st410_amr"

#' Data: Plasmid replicons identified in global ST167 isolates
#'
#' Identified plasmid replcions using ARIBA and the PlasmidFinder database.
#' See manuscript for details
#'
#' @format A tidy data frame with 1838 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{ref_seq}{Identified Plasmid replicon}
#' }
"btESBL_ecoli_st167_plasmids"


#' Data: Plasmid replicons identified in global ST410 isolates
#'
#' Identified plasmid replcions using ARIBA and the PlasmidFinder database.
#' See manuscript for details
#'
#' @format A tidy data frame with 3308 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{ref_seq}{Identified Plasmid replicon}
#' }
"btESBL_ecoli_st410_plasmids"

#' Data: Covariate data to pass to Stan longitudinal carriage model in list format
#'
#' See vignette for details. Each row encodes two ESBL observations a time
#' t apart and the covariate values in this time period
#'
#' These data are set up to fit a model where the effect of antimicrobial
#' exposure has a time-varying exponential decay effect, and hospitalisation
#' has a time-varying stepwise constant effect
#'
#' @format A named list of seven items.
#' \describe{
#'   \item{N}{Number of rows of data}
#'   \item{t}{vector of length n, Time between ESBL observations for that row}
#'   \item{n_covs}{Number of covariates to include in model}
#'   \item{covs_type}{define types of covariates}
#'   \item{start_state}{Starting state for row ESBL present absent}
#'   \item{end_state}{End state for that time period}
#'   \item{covs_mat}{Covariate values ofr each included covariate}
#' }
"btESBL_stanmodeldata"

#' Data: Covariate data to pass to Stan longitudinal carriage model in dataframe format
#'
#' See vignette for details. Each row encodes two ESBL observations a time
#' t apart and the covariate values in this time period
#'
#' These data are set up to fit a model where the effect of antimicrobial
#' exposure has a time-varying exponential decay effect, and hospitalisation
#' has a time-varying stepwise constant effect
#'
#' @format The data from btESBL_stanmodeldata in dataframe format
"btESBL_modeldata"



#' Data: metadata of global E. coli ST131 isolates
#'
#' Metadata of a collection of 926 E. coli ST131 genomes
#' including 862 from a study of global E. coli ST131
#' (Mcnally et at 2019 https://doi.org/10.1128/mBio.00644-19)
#' and the remainder sequenced as part of this study
#'
#'@format A dataframe with 926 rows and 4 variables
#'\describe{
#'   \item{Accession_number}{Unique sample-sequencing run ID or accession}
#'   \item{Year}{Year of sample collection}
#'   \item{Country}{Country of sample collection}
#'   \item{clade}{ST131 clade as defined in Mcnally et al}
#' }
"btESBL_ecoli_st131_metadata"


#' Data: AMR genes identified in global ST131 isolates
#'
#' Identified AMR genes using ARIBA and the SRST database.
#' See manuscript for details
#'
#' @format A tidy data frame with 5494 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{gene}{Identified AMR gene sequence}
#' }
"btESBL_ecoli_st131_amr"


#' Data: Plasmid replicons identified in global ST131 isolates
#'
#' Identified plasmid replcions using ARIBA and the PlasmidFinder database.
#' See manuscript for details
#'
#' @format A tidy data frame with 1838 rows and 2 variables:
#' \describe{
#'   \item{name}{Unique sample-sequencing run ID or accession}
#'   \item{ref_seq}{Identified Plasmid replicon}
#' }
"btESBL_ecoli_st131_plasmids"


#' Data: Contig cluster assignment in sensitivity analyses
#'
#' Contig cluster assignment varying cd-hit sequence identity from
#' 0.95 - 1.00 and length cutoff from 0-0.8
#' Care! The cluster ids are not unqiue! Only the
#' cluster id-len_diff_cutoff-ident_cutoff combination is unique. See
#' the code for how this is managed in practice.
#' See manuscript for further details
#'
#' @format A tidy data frame with 6426 rows and 8 variables:
#' \describe{
#'   \item{cluster}{Cluster id from cd-hit}
#'   \item{id}{Sequence id within cluster from cd-hit}
#'   \item{length}{Length of sequence, bases}
#'   \item{contig}{Contid unique identfier}
#'   \item{lane}{sample ID}
#'   \item{len_diff_cutoff}{Length difference cutoff (-s in cd-hit)}
#'   \item{ident_cutoff}{Sequence identity cutoff (-c in cd hit)}
#' }
"btESBL_contigclusters_sensax"


#' Data: Contig cluster multiple sequence alignment .paf files
#'
#' Multiple sequence alignments of ten largest contig clusters, generated
#' using minimap2, mapping all cluster sequences to the longest sequence (i.e.
#' the cluster representative sequenc from cd hit).
#'
#' @format A named list of ten data frames, where name is cluster id. Columns are:
#' \describe{
#'   \item{qname}{Query sequence name}
#'   \item{qlen}{Query sequence length, bases}
#'   \item{qstart}{Query start coordinate (0-based)}
#'   \item{qend}{Query end coordinate (0-based)}
#'   \item{strand}{‘+’ if query/target on the same strand; ‘-’ if opposite}
#'   \item{tname}{Target sequence name}
#'   \item{tlen}{Target sequence length, bases}
#'   \item{tstart}{Target start coordinate on orignal strand}
#'   \item{tend}{Target end coordinate on orignal strand}
#'   \item{nmatch}{Number of matching bases in the mapping}
#'   \item{alen}{Number of bases, including gaps, in the mapping}
#'   \item{mapq}{Mapping quality (0-255, with 255 if missing)}
#' }
"btESBL_contigclusters_msa_paf_files"


#' Data: Contig cluster multiple sequence alignments as PopGenome GENOME objects
#'
#' Multiple sequence alignments of ten largest contig clusters, generated
#' using minimap2, mapping all cluster sequences to the longest sequence (i.e.
#' the cluster representative sequenc from cd hit). Sequence alignments have
#' been generated as .fasta files from the minimap2 .SAM files, then the
#' fastas loaded using the PopGenome::readData function.
#'
#' @format A named list of ten data GENOME objects, where name is cluster id.
#'
"btESBL_contigclusters_msa_alignments"

#' Data: Antimicrobial sensitivity testing for a suvet of isolates
#'
#' Results from disc antimicrobial sensitivity testing (AST) according
#' to BSAC guidelines for six antimicrobials. Only a subset of isolates
#' underwent AST: 442 E. coli and 167 K. pnemoniae sequence complex (KpSC).
#'
#'@format A dataframe with 609 rows and 8 variables
#'\describe{
#'   \item{supplier_name}{Unique sample ID}
#'   \item{Organism}{E. coli or KpSC}
#'   \item{amikacin}{AST result to amikacin}
#'   \item{chloramphenicol}{AST result to chloramphenicol}
#'   \item{ciprofloxacin}{AST result to ciprofloxacin}
#'   \item{cotrimoxazole}{AST result to cotrimoxazole}
#'   \item{gentamicin}{AST result to gentamicin}
#'   \item{meropenem}{AST result to meropenem}
#' }
"btESBL_AST"
