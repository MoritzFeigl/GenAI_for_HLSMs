#!/bin/sh
# GenAI_para mHM optimization mHM runs
# define working directory
PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/

# tf update
cd $PROJECTPATH/mhm/deps/mpr
/home/fs71468/mfeigl/miniconda3/envs/GenAI_para_mhm_py3/bin/python3 -m src_python.pre_proc.update_tfs_in_fortran_source -c $PROJECTPATH/Training/6335125/config/mpr.nml -p <DATA_DIR>05_mhm/master_mhm_config/mpr_global_parameter_mhm.nml --clean

# Load modules
source ~/env/spackenv
export LD_LIBRARY_PATH=$LIBRARY_PATH
export FC=ifort
# compile mhm
now=$(date +"%T")
echo -e "Compiling time start: $now"
cd $PROJECTPATH/mhm/build
make -j 8  > /dev/null 2>&1
now=$(date +"%T")
echo -e "Compiling time end: $now"

# copy mhm into basin folders
cd $PROJECTPATH/Training
for d in */
do
  ln -sf $PROJECTPATH/mhm/build/mhm-0.5.9 $PROJECTPATH/Training/${d}config/mhm-0.5.9
  ln -sf $PROJECTPATH/mhm/build/mhm $PROJECTPATH/Training/${d}config/mhm
done
