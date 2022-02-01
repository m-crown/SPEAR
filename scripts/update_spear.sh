#|/usr/bin/env bash -eu
# Matt Bashton 2021
# Updates SPEAR
# Assumes spear conda env is already active is script is envoked by spear itself

# We need a command to run
[ $# -ne 1 ] && { echo -en "Error, usage: $(basename $0) [spear|all-data|all]\n\n" ; exit 1; }

BASE=$PWD
COMMAND=${1}

if [[ $COMMAND == "spear" ]]; then
    echo "Updating SPEAR repo"
    cd ${CONDA_PREFIX}/repo
    git pull
    cd ${BASE}

    # Install SPEAR files to data/ and bin/
    cp ${CONDA_PREFIX}/repo/data/* ${CONDA_PREFIX}/data
    cp ${CONDA_PREFIX}/repo/scripts/* ${CONDA_PREFIX}/bin
    cp ${CONDA_PREFIX}/repo/images/* ${CONDA_PREFIX}/images
    chmod +x ${CONDA_PREFIX}/bin/*.sh
    chmod +x ${CONDA_PREFIX}/bin/*.py
    cp ${CONDA_PREFIX}/repo/spear ${CONDA_PREFIX}/bin
    chmod +x ${CONDA_PREFIX}/bin/spear


elif [[ $COMMAND == "all-data" ]]; then
    echo "Updating all data sets"
    cd ${CONDA_PREFIX}/repo
    git -C data pull
    cd ${BASE}
     # Install SPEAR files to data/
    cp ${CONDA_PREFIX}/repo/data/* ${CONDA_PREFIX}/data

    # Download data from GitHub
    wget -q https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf -O ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf
    bgzip -f ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf.gz
    tabix -p vcf ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf.gz
    wget -q https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS/master/results/single_mut_effects/single_mut_effects.csv -O ${CONDA_PREFIX}/data/single_mut_effects.csv
    wget -q https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv -O ${CONDA_PREFIX}/data/escape_calculator_data.csv
    wget -q https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/bindingcalculator.py -O ${CONDA_PREFIX}/bin/bindingcalculator.py
    wget -q https://raw.githubusercontent.com/nataliateruel/data_Spike/main/vibentropy_occupancy_dmsdata.csv -O ${CONDA_PREFIX}/data/vibentropy_occupancy_dmsdata.csv


elif [[ $COMMAND == "all" ]]; then
    echo "Updating SPEAR, data sets and dependencies"
    cd ${CONDA_PREFIX}/repo
    git pull
    cd ${BASE}

    # Update the env
    CONDA_PATH=$(conda list -n base | awk 'NR == 1' | grep -oP '/\S+[^:]')
    source ${CONDA_PATH}/etc/profile.d/conda.sh
    conda config --set channel_priority strict
    # Update and remove old dependencies no longer in yml
    conda env update -n spear -q -f ${CONDA_PREFIX}/repo/environment.yml --prune

    # Download SnpEff and SnpSift
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O ${CONDA_PREFIX}/snpEff_latest_core.zip
    unzip -q -o ${CONDA_PREFIX}/snpEff_latest_core.zip -d ${CONDA_PREFIX}
    rm ${CONDA_PREFIX}/snpEff_latest_core.zip
    ${CONDA_PREFIX}/snpEff/data
    unzip -q -o ${CONDA_PREFIX}/repo/data/NC_045512.2.zip -d ${CONDA_PREFIX}/snpEff/data
    
    # Detect OS
    if [[ ${OSTYPE} == "darwin"* ]]; then
	OS="macOS"
    elif [[ ${OSTYPE} == "linux"* ]]; then
	OS="linux"
    else
	echo "Unsupported OS: ${OSTYPE}"
	exit 1
    fi

    # OS Dependent download
    if [[ ${OS} = "linux" ]]; then
	wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToVcf -O ${CONDA_PREFIX}/bin/faToVcf
	chmod +x ${CONDA_PREFIX}/bin/faToVcf
	wget -q https://github.com/brentp/vcfanno/releases/download/v0.3.3/vcfanno_linux64 -O ${CONDA_PREFIX}/bin/vcfanno
	chmod +x ${CONDA_PREFIX}/bin/vcfanno
    elif [[ ${OS} = "macOS" ]]; then
	wget -q http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/faToVcf -O ${CONDA_PREFIX}/bin/faToVcf
	chmod +x ${CONDA_PREFIX}/bin/faToVcf
	wget -q https://github.com/brentp/vcfanno/releases/download/v0.3.3/vcfanno_osx -O ${CONDA_PREFIX}/bin/vcfanno
	chmod +x ${CONDA_PREFIX}/bin/vcfanno
    fi

    # Download data from GitHub
    wget -q https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf -O ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf
    bgzip -f ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf.gz
    tabix -p vcf ${CONDA_PREFIX}/data/problematic_sites_sarsCov2.vcf.gz
    wget -q https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS/master/results/single_mut_effects/single_mut_effects.csv -O ${CONDA_PREFIX}/data/single_mut_effects.csv
    wget -q https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv -O ${CONDA_PREFIX}/data/escape_calculator_data.csv
    wget -q https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/bindingcalculator.py -O ${CONDA_PREFIX}/bin/bindingcalculator.py
    wget -q https://raw.githubusercontent.com/nataliateruel/data_Spike/main/vibentropy_occupancy_dmsdata.csv -O ${CONDA_PREFIX}/data/vibentropy_occupancy_dmsdata.csv

    # Install SPEAR files to data/ and bin/
    cp ${CONDA_PREFIX}/repo/data/* ${CONDA_PREFIX}/data
    cp ${CONDA_PREFIX}/repo/scripts/* ${CONDA_PREFIX}/bin
    chmod +x ${CONDA_PREFIX}/bin/*.sh
    chmod +x ${CONDA_PREFIX}/bin/*.py
    cp ${CONDA_PREFIX}/repo/spear ${CONDA_PREFIX}/bin
    chmod +x ${CONDA_PREFIX}/bin/spear

    # Test snpEff is working
    java -Xmx1g -jar $CONDA_PREFIX/snpEff/snpEff.jar -hgvs1LetterAa -nodownload -no SPLICE_SITE_ACCEPTOR -no SPLICE_SITE_DONOR -no SPLICE_SITE_REGION -no SPLICE_SITE_BRANCH -no SPLICE_SITE_BRANCH_U12 -noLog -noLof -no-intron -noMotif -noStats -no-downstream -no-upstream -no-utr NC_045512.2 $CONDA_PREFIX/data/install_vcf.vcf > /dev/null

fi
