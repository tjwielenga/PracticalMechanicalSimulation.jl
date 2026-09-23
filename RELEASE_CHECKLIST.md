# Public Release Checklist

This checklist records the remaining work for the first public release. It is
not a substitute for the living [Project Summary](PROJECT_SUMMARY.md).

## Repository

- [x] Choose and add the repository license (MIT).
- [x] Create the public repository from a clean release snapshot and configure
      its Git remote. Do not publish the private development history, which
      contains the excluded research archive in earlier commits.
- [x] Add the final repository URL to the README and citation metadata.
- [x] Confirm that no credentials, private paths, generated results, or Word
      lock files are tracked.
- [x] Keep the author's local research archive outside the public repository.

## Software

- [x] Run the complete supported Julia test suite from a clean detached
      release snapshot.
- [x] Run the SimpView tests and production build from that clean snapshot.
- [x] Confirm the command-line quick start on macOS.
- [x] Confirm the command-line quick start on Linux through public CI.
- [x] Regenerate representative planar and spatial results with release code.
- [x] Keep SimpView as the sole supported viewer and remove superseded
      viewer dependencies.

## Documentation

- [x] Check all local Markdown links.
- [x] Read the README as a new user and repeat its installation procedure from
      a fresh clone of the public `v0.1.0` tag.
- [x] Confirm that planar and spatial model references match the implemented
      TOML fields.
- [x] Confirm that public limitations are stated in the README and Project
      Summary.
- [x] Add the public repository URL and release identifier to the citation
      file.

## Release

- [x] Set the release version in `Project.toml` and SimpView metadata.
- [x] Record release notes and known limitations.
- [x] Tag the tested commit as `v0.1.0`.
- [x] Publish the repository and verify the macOS, Linux, and SimpView
      continuous-integration checks.
- [ ] Archive or link the exact supporting models and results used by the
      methods paper when they are ready for publication.

## Latest local validation

September 16, 2026:

- clean-snapshot supported Julia suite: 2,345 of 2,345 checks passed;
- SimpView: 19 of 19 checks passed and the production bundle built;
- npm production dependency audit: no known vulnerabilities reported;
- public CI: SimpView and Julia 1.12 on macOS and Ubuntu passed;
- model-reference audit: every registered planar and spatial element type is
  covered, 70 standalone distributed models load, and the two remaining
  continuation models correctly request their prerequisite static results;
- planar quick start: complete 81-sample kinematic result;
- spatial quick start: complete 301-sample dynamic result; and
- Large Van: six-second, 361-frame dynamic comparison at 30 m/s with normal
  and raised center of mass and roof-ground contacts.
