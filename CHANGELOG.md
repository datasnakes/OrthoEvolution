# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive test suites for `utilities`, `Manager`, and `Tools` modules
  - Added 20+ test methods for `BlastUtils`, `GenbankUtils`, `OrthologUtils`, `ManagerUtils`, `PackageVersion`, and `FullUtilities` classes
  - Added 16 test methods for `Management`, `RepoManagement`, `UserManagement`, `WebsiteManagement`, `ProjectManagement`, and `BaseDatabaseManagement` classes
  - Added 5 test methods for `LogIt`, `PyBasher`, `NcbiFTPClient`, and `Multiprocess` classes
  - Total of ~705 lines of new test code across three test files
- Docstrings added to Orthologs module

### Fixed
- Fixed typo in `FullUtilities.__init__()`: `Ortho0logUtils` → `OrthologUtils` (resolves `AttributeError`)
- Fixed missing docstrings across codebase
- Fixed various test issues and code coverage gaps

### Changed
- Updated Python version requirements and examples
- Refactored code structure (#177)
- Moved templates to `.github` directory
- Updated Python tasks configuration

## [1.0.0b1] - Development

### Added
- Test infrastructure improvements and codecov integration
- GitHub Actions CI/CD pipeline replacing Travis CI
- pytest configuration and test suite expansion
- Test GPCR dataset for testing
- PyBasher commands for bash operations
- Tarfile member sanitization for security (#173)

### Changed
- Removed Travis CI in favor of GitHub Actions
- Updated Python version support (3.9, 3.10, 3.11)
- Upgraded Biopython version
- Bumped dependencies: tqdm (4.25.0 → 4.66.3), jinja2, werkzeug, sqlalchemy

### Fixed
- Fixed tarfile member sanitization in extractall() (#173)
- Fixed various test failures and code coverage issues
- Fixed requirements.txt security vulnerabilities (#179, #175, #176, #170, #169, #168)
- Fixed Cookies tests and Oven/CookBook text issues
- Fixed log level color formatting
- Fixed Flask version requirements and psutil requirement

## [0.9.0a2] - Previous Release

### Added
- Initial test suite
- Basic CI/CD setup

### Changed
- Package structure and organization

---

## Version History Notes

- **1.0.0b1**: Current development version (beta)
- **0.9.0a2**: Previous tagged release (alpha)
