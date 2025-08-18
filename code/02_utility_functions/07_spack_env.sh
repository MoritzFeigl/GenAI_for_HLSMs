#!/bin/sh
#
# Generate module load commands in ~/env/spackenv
mkdir $HOME/env
cat <<EOF | /bin/sh >$HOME/env/spackenv
FIND='spack module tcl loads --dependencies'
export MODULES_AUTO_HANDLING=1
\$FIND intel-oneapi-compilers/56q5kyk
\$FIND intel-oneapi-mpi/ddvuxhy
\$FIND parallel-netcdf/vtmjge7
\$FIND netcdf-c/y2uuoro
\$FIND hdf5/2yy5m3y
\$FIND netcdf-fortran/jrv4ciq
\$FIND cmake/cfhkzu5
\$FIND r/d3w6clx
EOF

# old version
#\$FIND intel-oneapi-compilers@2022.0.2
#\$FIND intel-oneapi-mpi@2021.4.0%intel@2021.5.0
#\$FIND /rsjx4ta
#\$FIND /6xyj236
#\$FIND netcdf-fortran
#\$FIND parallel-netcdf@1.12.2
#\$FIND cmake@3.21.4%intel@2021.5.0
#\$FIND r
# newrer old version
#\$FIND intel-oneapi-compilers@2022.0.2
#\$FIND intel-oneapi-mpi@2021.4.0%intel@2021.5.0
#\$FIND hdf5/tty2bao
#\$FIND netcdf-c/spzlhyr
#\$FIND parallel-netcdf/kuvf5h6
#\$FIND netcdf-fortran/um5yjit
#\$FIND cmake@3.21.4%intel@2021.5.0
#\$FIND r
#EOF
