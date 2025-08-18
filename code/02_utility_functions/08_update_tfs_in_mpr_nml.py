import f90nml
import pandas as pd
import argparse as ap
import os
from multiprocessing import Pool

if __name__ == "__main__":
    # script to update TF definitions in mpr.nml
    parser = ap.ArgumentParser()
    parser.add_argument('-f', metavar='--file',
                        required=True, help='Path to basin directory, where each basin has a config/mpr.nml')
    parser.add_argument('-t', metavar='--tfs',
                        required=True, help='Path to tf csv file.')
    parser.add_argument('-s', metavar='--save',
                        required=True, help='Save parameters as ncdf?.')

    # get paths
    args = parser.parse_args()

    # get a mpr file from f directory
    all_basins = os.listdir(args.f)
    first_basin = all_basins[1]

    parser = f90nml.Parser()
    parser.global_start_index = 1
    nml = parser.read(f"{args.f}/{first_basin}/config/mpr.nml")
    tfs = pd.read_csv(args.t)

    # change all parameter TFs and data_array lists for all basins
for i in range(tfs.shape[0]):
    print(i)

    para_name = tfs["parameters"][i]
    para_tf = tfs["TFs"][i]
    par_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == para_name][0]
    nml['data_arrays']['transfer_func'][par_ind] = para_tf
    nml['data_arrays']["from_data_arrays"][par_ind] = tfs["inputs"][i].split(", ")

    if args.s:
        important_paras = ["KSat", "KSat_till", "KSat_notill", "vGenu_n_till", "vGenu_n_notill", "ThetaS_till",
                           "ThetaS_notill", "fRoots_temp", "L1_Max_Canopy_Intercept", "FieldCap_till", "FieldCap_notill"]
        for para_name in important_paras:
            par_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == para_name][0]
            nml['data_arrays']['to_file'][par_ind] = True

    # overwrite mpr.nml
    def task(basin):
      nml_basin = parser.read(f"{args.f}/{basin}/config/mpr.nml")
      nml['data_arrays']["from_file"] = nml_basin['data_arrays']["from_file"]
      nml.write(f"{args.f}/{basin}/config/mpr.nml", force=True)

    with Pool() as pool:
        # call the function for each item in parallel
        pool.map(task, all_basins)

    print("Updated all parameter TFs.\n")
