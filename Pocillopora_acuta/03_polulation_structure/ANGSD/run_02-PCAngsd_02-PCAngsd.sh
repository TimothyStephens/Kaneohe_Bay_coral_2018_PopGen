#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py38; set -eu

NCPUS=48
BEAGLE="PCAngsd.angsd.beagle.gz"

#### Start Script
# See http://www.popgen.dk/software/index.php/PCAngsdTutorial

# Estimating Individual Allele Frequencies
run_cmd "pcangsd --threads ${NCPUS} --beagle ${BEAGLE} --out ${BEAGLE}.IndAlleleFreq"

# Without Estimating Individual Allele Frequencies
run_cmd "pcangsd --threads ${NCPUS} --beagle ${BEAGLE} --out ${BEAGLE}.WithOutIndAlleleFreq --iter 0"

# Admixture based on two PC
run_cmd "pcangsd --threads ${NCPUS} --beagle ${BEAGLE} --out ${BEAGLE}.Admixture --admix --admix_alpha 50"


