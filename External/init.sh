#!/bin/sh

PYTHIA_VERSION=6428

PYTHIA_DIR=./
PYTHIA_FILE=pythia${PYTHIA_VERSION}.f
if [ ! -f ${PYTHIA_DIR}${PYTHIA_FILE} ]; then
  curl "https://pythia.org/download/pythia6/${PYTHIA_FILE}" > ${PYTHIA_DIR}${PYTHIA_FILE}
fi

HERWIG_VERSION=6521
HERWIG_DIR=./
declare -a HERWIG_FILES=("herwig${HERWIG_VERSION}.f" "herwig${HERWIG_VERSION}.inc" "HERWIG65.INC")
for HERWIG_FILE in "${HERWIG_FILES[@]}"; do
  if [ ! -f ${HERWIG_DIR}${HERWIG_FILE} ]; then
    curl "https://www.hep.phy.cam.ac.uk/theory/webber/Herwig/${HERWIG_FILE}" > ${HERWIG_DIR}${HERWIG_FILE}
  fi
done

if [[ $(hostname -f) =~ ^lxplus[0-9]+.cern.ch ]]; then
  echo ">>> on LXPLUS"
  curl "https://raw.githubusercontent.com/BristolTopGroup/AnalysisSoftware/master/FindROOT.cmake" > cmake/FindROOT.cmake
fi
