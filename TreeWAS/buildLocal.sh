#!/bin/bash

PKG_NAME=$( cat DESCRIPTION | grep "Package:" | sed -e "s/Package: //g")
Rscript -e "roxygen2::roxygenize(\"../${PKG_NAME}\")"

# Building package
R CMD build  .

PKG_VERSION=$( cat DESCRIPTION | grep Version | sed -e "s/Version: //g")
PKG_TARBALL=${PKG_NAME}_${PKG_VERSION}.tar.gz
RCHECK_DIR=${PKG_NAME}.Rcheck

# Checking package
R CMD check "${PKG_TARBALL}" --as-cran; CHECK_RET=r-devel

# R CMD check fail logs
for name in $(find "${RCHECK_DIR}" -type f -name "*fail");do echo ">>> Filename: ${name} <<<";cat ${name};done

# Uncomment the following line to show R CMD check log logs
#for name in $(find "${RCHECK_DIR}" -type f -name "*log");do echo ">>> Filename: ${name} <<<";cat ${name};done

# Uncomment the following line to show R CMD check out logs
#for name in $(find "${RCHECK_DIR}" -type f -name "*out");do echo ">>> Filename: ${name} <<<";cat ${name};done

