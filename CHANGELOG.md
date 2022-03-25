# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2022-03-25
### Added
- Longitudinal carriage model ceftriaxone vs all other abx
- btESBL_AST: phenotypic AST data and plots in carriage paper
- btESBL_model2simulations_2: further simulations from posterior
- Contig cluster sensitivity analysis and mapping of ID/AMR genes/plasmid replicon 

### Fixed
- bTESBL_plasmidreplicons didn't include kleb in error - fixed

### Changed
- Correct typos throughout

## [1.1.4] - 2021-12-13
### Added
- Table of Kleb isolates showing ESBL status and which study they came from
- Within-participant Kleb ST correlation plots
- Sensitivity analysis plots
- Add VFDB virulence screening for Kleb

## [1.1.3] - 2021-11-12
### Changed
- Tweak E. coli plots

## [1.1.2] - 2021-11-12
### Added
- Add documentation of E. coli ST 131 data

## [1.1.1] -  2021-11-12
### Changed
- Switch E coli Malawi carriage tree to non ASC 
- Switch E coli global tree to non ASC, using full Horesh collection
- Use ClermonTyping v20.3 instead of isPcr for phylogroups
- switch E coli ST410 and ST167 trees to non-ASC
- add E coli ST131 tree and metadata
### Fixed 
- Updated kleb Malawi carriage tree to non ASC - 1.1.0 didn't have correct ones


## [1.1.0] - 2021-11-01
### Changed
- Switch Klebsiella Malawi carriage tree to be without ASC
- Switch Kleb global tree to be without ASC 
- Switch Kleb Malawi all isolate tree to be without ASC
- Switch to semantic versioning

## [1.0.0.3] - 2021-10-12
### Fixed
- remove second geom_treescale on global klebs tree plot
- change heatmap on malawi kleb tree for better pdf renderning

## [1.0.0.2] - 2021-10-11
### Fixed
- Add rstan dependency

## [1.0.0.1] - 2021-10-11
### Fixed
- Add ggtree bioconductor dependency with Remotes:
- Add missing quote to readme install_github call



