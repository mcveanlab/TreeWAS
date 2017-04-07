#!/bin/bash
# Building package
R CMD build  .

TREEWAS_VERSION=$( cat DESCRIPTION | grep Version | sed -e "s/Version: //g")
PKG_TARBALL=TreeWAS_${TREEWAS_VERSION}.tar.gz

# Checking package
R CMD check "${PKG_TARBALL}" --as-cran; CHECK_RET=r-devel

# R CMD check fail logs
$ for name in $(find "${RCHECK_DIR}" -type f -name "*fail");do echo ">>> Filename: ${name} <<<";cat ${name};done

# R CMD check log logs
$ for name in $(find "${RCHECK_DIR}" -type f -name "*log");do echo ">>> Filename: ${name} <<<";cat ${name};done

# R CMD check out logs
$ for name in $(find "${RCHECK_DIR}" -type f -name "*out");do echo ">>> Filename: ${name} <<<";cat ${name};done

