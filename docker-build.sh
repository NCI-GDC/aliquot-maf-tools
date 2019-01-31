#!/bin/bash
# This script does the following:
# - Overwrites .deps dir TODO: Only over write required repos only if not on hash
# - Parses requirements.txt and clones all private repos in dependencies
# - Creates docker.requirements.txt with remaining modules

# Overwriting .deps dir
D_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.deps"
[[ -d ${D_DIR} ]] && rm -rf ${D_DIR}/* || mkdir -p ${D_DIR}

pushd () {
    command pushd "$@" > /dev/null
}

popd () {
    command popd "$@" > /dev/null
}

# Parses requirements.txt and clones all private repos in dependencies to a specific hash
for i in $(grep '^-e\|\.git' requirements.txt | sed -n 's/.*git+\([^ ]*\)#egg.*/\1/p')
do
  echo "Git Repo:${i%@*} Git Hash:${i##*@}"
  pushd ${D_DIR}
  repo=${i%@*}
  git clone ${repo}
  pushd `__=${repo%.git}&&echo ${__##*/}`
  echo $PWD
  git reset --hard ${i##*@}
  popd
  popd
done

# Creating requirements file with remaining modules
grep -v "^-e\|\.git" requirements.txt > docker-requirements.txt

# tag
quay="quay.io/ncigdc/aliquot-maf-tools"
version="latest"
imagetag="${quay}:${version}"

echo "Building tag: $imagetag"
docker_build -t $imagetag .

# Cleanup
rm docker-requirements.txt
if [[ -d ${D_DIR} ]]; then rm -rf ${D_DIR}; fi
