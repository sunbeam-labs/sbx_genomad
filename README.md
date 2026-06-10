<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_genomad

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_genomad/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_genomad/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_genomad)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_genomad/)
<!-- badges: end -->

## Introduction

sbx_genomad is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for identifying viruses in samples with [geNomad](https://portal.nersc.gov/genomad/index.html). This pipeline uses [MEGAHIT](https://github.com/voutcn/megahit) for assembly of contigs and then processes assemblies with geNomad.

N.B. This extension requires also having sbx_assembly installed.

### Installation

```
sunbeam extend https://github.com/sunbeam-labs/sbx_assembly.git
sunbeam extend https://github.com/sunbeam-labs/sbx_genomad.git
```

### Database

sbx_genomad expects the genomad reference database to be available locally. Download the database following the official instructions, for example:

```
genomad download-database .
```

Update the `genomad_db` entry in your Sunbeam configuration to point at the resulting directory.

### Running

Run with sunbeam on the target `all_genomad`:

```
sunbeam run --profile /path/to/project/ --skip decontam all_genomad
```

### Options for config.yml

  - genomad_db: path to geNomad db (default: "/mnt/isilon/hvp/dbs/genomad_db") (NOTE: this should be a directory)
