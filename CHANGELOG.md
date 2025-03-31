# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.4.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs.
- Path to STAR index.
- Module to run STAR.
- Flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata.

### Changed
- Changes to medata.

## [2.3.0] - 2023-06-13
### Change
- Move assembly-specific settings into the wdl.

## [2.2.0] - 2023-02-23
### Change
- Update STAR to 2.7.10b.
- Change HardClip to Softclip for Chimeric Reads.

## [2.0.4] - 2022-06-24
### Change
- Update calculate.sh for the junction file.
- Updated README.md

## [2.1.1] - 2021-05-31
### Change
- Migrate to Vidarr.

## [2.1] - 2021-01-14
### Changed
- Update STAR to 2.7.6a.

### Added
- Add chimOutType option to write chimeric reads WithinBAM.

## [2.0.2] - 2020-05-04
### Added
- Make --chimOutJunctionFormat parameter optional for versions of star older than 2.6.1a.

## [2.0.1] - 2020-01-13
### Added
- Hotfix to expose --chimOutJunctionFormat parameter needed for starFusion.

## [2.0] - 2019-09-29
### Changed
- Converting to cromwell wdl with changes to provision junction file and take multiple fastq pairs.

## [1.1.1] - 2017-07-12
### Fixed
- Fixed a bug with adapter sequences paramtrizing.

## [1.1] - 2016-10-03
### Added
- Added new parameter produce_transcriptome_bam (true by default).

## [1.0] - 2016-09-08
### Added
- Initial import of STAR code.
