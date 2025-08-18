import f90nml
import pandas as pd
import argparse as ap
import os
mpr = "/mnt/Data/Dropbox/Projekte/GenAI_para_for_HLSMs/data/05_mhm/master_mhm_config/mpr_original.nml"


# get a mpr file from f directory
parser = f90nml.Parser()
parser.global_start_index = 1
nml = parser.read(mpr)

# fix limit length
add_lim = len(nml['data_arrays']["name"]) - len(nml['data_arrays']["limits"])
for add in range(add_lim):
    nml['data_arrays']["limits"] = nml['data_arrays']["limits"] + [[None, None]]
add_tf = len(nml['data_arrays']["name"]) - len(nml['data_arrays']["transfer_func"])
for add in range(add_tf):
    print(add)
    nml['data_arrays']["transfer_func"] = nml['data_arrays']["transfer_func"] + [None]

# change lai source and move forward
lai_id = [i for i, x in enumerate(nml['data_arrays']['name']) if x == "lai_class"][0]
nml['data_arrays']['from_file'][lai_id] = '<DATA_DIR>06_training_validation_data/<SPLIT>/<BASIN>/static/mpr/lai.nc'
new_lai_id = 1
arrays = ['name', 'from_file', 'to_file', 'limits', 'from_data_arrays', 'target_coord_names', 'upscale_ops', 'transfer_func']
for array in arrays:
    values = nml['data_arrays'][array]
    values = values[:new_lai_id] + [values[lai_id]] + values[new_lai_id:lai_id] + values[(lai_id+1):]
    nml['data_arrays'][array] = values

# add lai_agg
add_var = ["lai_agg"]
start_ind = 1
nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][start_ind:]
nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
                                  ['<DATA_DIR>06_training_validation_data/<SPLIT>/<BASIN>/static/mpr/lai_agg.nc'] \
                                  + nml['data_arrays']['from_file'][start_ind:]
nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
                                  [False] + nml['data_arrays']['to_file'][start_ind:]
nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
[[None, None]] + nml['data_arrays']['limits'][start_ind:]
nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
                               [[None, None, None]] +\
                                         nml['data_arrays']['from_data_arrays'][start_ind:]
nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
                                         [[None, None]] + \
                                           nml['data_arrays']['target_coord_names'][start_ind:]
upscale_template = [None, None, None]
nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
                                    [upscale_template] + \
                                    nml['data_arrays']['upscale_ops'][start_ind:]
nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
                                    [None] + nml['data_arrays']['transfer_func'][start_ind:]


# move mat, mat_range, map forward
for var in ["mat", "mat_range", "map"]:
    var_id = [i for i, x in enumerate(nml['data_arrays']['name']) if x == var][0]
    new_id = 2
    arrays = ['name', 'from_file', 'to_file', 'from_data_arrays', 'target_coord_names', 'upscale_ops', 'transfer_func']
    for array in arrays:
        values = nml['data_arrays'][array]
        if var_id == len(values)-1:
            values = values[:new_id] + [values[var_id]] + values[new_id:var_id]
        else:
            values = values[:new_id] + [values[var_id]] + values[new_id:var_id] + values[(var_id + 1):]
        nml['data_arrays'][array] = values

    values = nml['data_arrays']["limits"]
    if var_id == len(values) - 1:
        values = values[:new_id] + [values[var_id]] + values[new_id:var_id]
    else:
        values = values[:new_id] + [values[var_id]] + values[new_id:var_id] + [values[(var_id + 1):]]
    nml['data_arrays']["limits"] = values

# check if length is always the same
#arrays = ['name', 'to_file', 'limits', 'from_data_arrays', 'target_coord_names', 'upscale_ops', 'transfer_func']
#for array in arrays:
#    print(array + str(len(nml['data_arrays'][array])))


# add _till and _notill variables
for addon in ['_notill', '_till']:
    for var in ['slope', 'dem', 'aspect', "mat", "mat_range", "map", "lai_agg"]:
        add_var = [var + addon]
        start_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == var][0]+1
        # names
        nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][
                                                                                         start_ind:]
        # from_file
        nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
                                          [None] + nml['data_arrays']['from_file'][start_ind:]
        # to file
        nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
                                        [False] + nml['data_arrays']['to_file'][start_ind:]
        # limits
        nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
                                       [[None, None]] + nml['data_arrays']['limits'][start_ind:]
        # from_data_arrays
        nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
                                                 [[var, None, None]] + \
                                                 nml['data_arrays']['from_data_arrays'][start_ind:]
        # target_coord_names
        nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
                                                   [['land_cover_period', 'lat', 'lon', 'horizon' + addon]] + \
                                                   nml['data_arrays']['target_coord_names'][start_ind:]
        # upscale_ops
        nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
                                            [['1.0', '1.0', '1.0', '1.0']] + \
                                            nml['data_arrays']['upscale_ops'][start_ind:]
        # transfer_func
        nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
                                              [None] + \
                                              nml['data_arrays']['transfer_func'][start_ind:]


add_var = ['slope_horizon', 'dem_horizon', 'aspect_horizon', 'mat_horizon', 'mat_range_horizon',
           'map_horizon', "lai_agg_horizon"]
start_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == "bd_eff_till"][0] + 1
nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][start_ind:]
nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
                                  [None, None, None, None, None, None, None] + nml['data_arrays']['from_file'][start_ind:]
nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
                                  [False, False, False, False, False, False, False] + nml['data_arrays']['to_file'][start_ind:]
nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
[[None, None], [None, None], [None, None], [None, None], [None, None], [None, None], [None, None]] + nml['data_arrays']['limits'][start_ind:]
nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
                               [['slope', None, None],
                                ['dem', None, None],
                                ['aspect', None, None],
                                 ['mat', None, None],
                                 ['mat_range', None, None],
                                 ['map', None, None],
                                ['lai_agg', None, None]] +\
                                         nml['data_arrays']['from_data_arrays'][start_ind:]
nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
                                         [['lat', 'lon', 'horizon'],
                                          ['lat', 'lon', 'horizon'],
                                          ['lat', 'lon', 'horizon'],
                                           ['lat', 'lon', 'horizon'],
                                           ['lat', 'lon', 'horizon'],
                                           ['lat', 'lon', 'horizon'],
                                          ['lat', 'lon', 'horizon']] + \
                                           nml['data_arrays']['target_coord_names'][start_ind:]
upscale_template = ['1.0', '1.0', '1.0']
nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
                                    [upscale_template, upscale_template, upscale_template,
                                     upscale_template, upscale_template, upscale_template, upscale_template] + \
                                    nml['data_arrays']['upscale_ops'][start_ind:]
nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
                                    [None, None, None, None, None, None, None] + nml['data_arrays']['transfer_func'][start_ind:]









# for fRoots
# use horizon mean for: sand, clay, bd
# add land_cover dim to all variables!
##### Add froots variables
for var in ['bd', 'sand', 'clay', "mat", "mat_range", "map", "lai_agg"]:
    add_var = [var + '_land_cover']
    start_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == var][0]+1
    # names
    nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][
                                                                                     start_ind:]
    # from_file
    nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
                                      [None] + nml['data_arrays']['from_file'][start_ind:]
    # to file
    nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
                                    [False] + nml['data_arrays']['to_file'][start_ind:]
    # limits
    nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
                                   [[None, None]] + nml['data_arrays']['limits'][start_ind:]
    # from_data_arrays
    nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
                                             [[var, None, None]] + \
                                             nml['data_arrays']['from_data_arrays'][start_ind:]
    # target_coord_names
    nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
                                               [['land_cover_period', 'lat', 'lon', 'horizon_out']] + \
                                               nml['data_arrays']['target_coord_names'][start_ind:]
    # upscale_ops
    nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
                                        [['1.0', '1.0', '1.0', '1.0']] + \
                                        nml['data_arrays']['upscale_ops'][start_ind:]
    # transfer_func
    nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
                                          [None] + \
                                          nml['data_arrays']['transfer_func'][start_ind:]

for var in ['slope', 'dem', 'aspect']:
    add_var = [var + '_land_cover']
    start_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == var + "_horizon"][0]+1
    # names
    nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][
                                                                                     start_ind:]
    # from_file
    nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
                                      [None] + nml['data_arrays']['from_file'][start_ind:]
    # to file
    nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
                                    [False] + nml['data_arrays']['to_file'][start_ind:]
    # limits
    nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
                                   [[None, None]] + nml['data_arrays']['limits'][start_ind:]
    # from_data_arrays
    nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
                                             [[var + "_horizon", None, None]] + \
                                             nml['data_arrays']['from_data_arrays'][start_ind:]
    # target_coord_names
    nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
                                               [['land_cover_period', 'lat', 'lon', 'horizon_out']] + \
                                               nml['data_arrays']['target_coord_names'][start_ind:]
    # upscale_ops
    nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
                                        [['1.0', '1.0', '1.0', '1.0']] + \
                                        nml['data_arrays']['upscale_ops'][start_ind:]
    # transfer_func
    nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
                                          [None] + \
                                          nml['data_arrays']['transfer_func'][start_ind:]


# add _horizon variables
#for var in ['slope', 'dem', 'aspect', "mat", "mat_range", "map", "lai_agg"]:
#    add_var = [var + '_horizon']
#    start_ind = [i for i, x in enumerate(nml['data_arrays']['name']) if x == var][0]+1
#    # names
#    nml['data_arrays']['name'] = nml['data_arrays']['name'][0:start_ind] + add_var + nml['data_arrays']['name'][
#                                                                                     start_ind:]
#    # from_file
#    nml['data_arrays']['from_file'] = nml['data_arrays']['from_file'][0:start_ind] + \
#                                      [None] + nml['data_arrays']['from_file'][start_ind:]
#    # to file
#    nml['data_arrays']['to_file'] = nml['data_arrays']['to_file'][0:start_ind] + \
#                                    [False] + nml['data_arrays']['to_file'][start_ind:]
#    # limits
#    nml['data_arrays']['limits'] = nml['data_arrays']['limits'][0:start_ind] + \
#                                   [[None, None]] + nml['data_arrays']['limits'][start_ind:]
#    # from_data_arrays
#    nml['data_arrays']['from_data_arrays'] = nml['data_arrays']['from_data_arrays'][0:start_ind] + \
#                                             [[var , None, None]] + \
#                                             nml['data_arrays']['from_data_arrays'][start_ind:]
#    # target_coord_names
#    nml['data_arrays']['target_coord_names'] = nml['data_arrays']['target_coord_names'][0:start_ind] + \
#                                               [['lat', 'lon', 'horizon']] + \
#                                               nml['data_arrays']['target_coord_names'][start_ind:]
#    # upscale_ops
#    nml['data_arrays']['upscale_ops'] = nml['data_arrays']['upscale_ops'][0:start_ind] + \
#                                        [['1.0', '1.0', '1.0']] + \
#                                        nml['data_arrays']['upscale_ops'][start_ind:]
#    # transfer_func
#    nml['data_arrays']['transfer_func'] = nml['data_arrays']['transfer_func'][0:start_ind] + \
#                                          [None] + \
#                                          nml['data_arrays']['transfer_func'][start_ind:]#
#



# write mpr.nml
new_filename = "/".join(mpr.split("/")[:-1] + ["mpr.nml"])
nml.write(new_filename, force=True)


[i for i, x in enumerate(nml['data_arrays']['name']) if x == "fRoots_temp"][0]
nml['data_arrays']['limits'][84]
[i for i, x in enumerate(nml['data_arrays']['name']) if x == "land_cover_horizon"][0]
nml['data_arrays']['target_coord_names'][83]

[i for i, x in enumerate(nml['data_arrays']['name']) if x == "fRoots_rescaled_total"][0]
nml['data_arrays']['transfer_func'][99]
