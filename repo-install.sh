#!/bin/bash
# This script installs python directories in the order
# of dependencies.

# Example: There are for python repos A to D in .deps dir
# if A Depends on A, B, C
# B -> C, D
# C -> D
# Then the script will install in the following order:
# D, C, B, A

pushd ${1:-.deps} > /dev/null #Getting inside .deps

# Getting list of repos inside .deps
__=($(ls -d */))
GLIST=(${__[*]%/})

# Function: install_module
# Arguments:
#  1. list of repos
#  2. repo to install
install_module () {
  declare -a list=("${!1}")
  local mod=${2}
  declare -a rem_list=(${list[@]/$mod})
  for r_mod in "${rem_list[@]}"
  do
    if grep -q ${r_mod} ${mod}/setup.py
    then
      install_module rem_list[@] ${r_mod}
    fi
  done
  if [[ " ${GLIST[@]} " =~ " ${mod} " ]]
  then
    pushd $mod > /dev/null
    pip install .
    popd > /dev/null
    GLIST=(${GLIST[@]/$mod})
  fi
  [[ ${#GLIST[@]} -ne 0 ]] && install_module GLIST[@] ${GLIST[0]}
}

# Calling function to install repos in order.
install_module GLIST[@] ${GLIST[0]}

popd > /dev/null #leaving .deps
