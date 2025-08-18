#!/bin/sh
# GenAI_para mHM optimization mHM runs

PROJECTPATH=/home/fs71468/mfeigl/GenAI_para_for_HLSMs/runs/

# update tfs
cd $PROJECTPATH/mhm/deps/mpr
/home/fs71468/mfeigl/miniconda3/envs/GenAI_para_mhm_py3/bin/python3 -m src_python.pre_proc.update_tfs_in_fortran_source -c $PROJECTPATH/Training/6335125/config/mpr.nml -p <DATA_DIR>05_mhm/master_mhm_config/mpr_global_parameter_mhm.nml --clean

# compile mhm
source ~/env/spackenv
export LD_LIBRARY_PATH=$LIBRARY_PATH
now=$(date +"%T")
echo -e "Compiling time start: $now"
cd $PROJECTPATH/mhm/
rm -r build
mkdir build
cd build
FC=ifort cmake -DCMAKE_WITH_OpenMP:STRING=OFF -DCMAKE_BUILD_TYPE=Debug -DNETCDF_LIBRARY=/gpfs/opt/sw/zen/spack-0.19.0/opt/spack/linux-almalinux8-zen3/intel-2021.7.1/netcdf-c-4.9.0-y2uuoroyrups2fhnjrqwzfyqz2tqvdpe/lib/libnetcdf.so -DNETCDF_F90_LIBRARY=/gpfs/opt/sw/zen/spack-0.19.0/opt/spack/linux-almalinux8-zen3/intel-2021.7.1/netcdf-fortran-4.6.0-jrv4ciqhozgatojzant24uxwnvstkw7s/lib/libnetcdff.so ..
make -j 8
now=$(date +"%T")
echo -e "Compiling time end: $now"
# copy mhm into basin folders

cd $PROJECTPATH/Training
for d in */
do
  ln -sf $PROJECTPATH/mhm/build/mhm-0.5.9 $PROJECTPATH/Training/${d}config/mhm-0.5.9
  ln -sf $PROJECTPATH/mhm/build/mhm $PROJECTPATH/Training/${d}config/mhm
done
