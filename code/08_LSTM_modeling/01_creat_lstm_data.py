import os
import xarray as xr
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
path = "C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data"

for split in ["Training", "Validation"]:
    split_dir = path + f"/06_training_validation_data/{split}/"
    split_basins = [name for name in os.listdir(split_dir) if os.path.isdir(os.path.join(split_dir, name))]

    # create basin-wise mean data frames

    # Data structure:
    lstm_data_dir = path + f"/09_lstm_data/"
    os.makedirs(lstm_data_dir, exist_ok=True)
    os.makedirs(lstm_data_dir + "/static", exist_ok=True)
    os.makedirs(lstm_data_dir + "/dynamic", exist_ok=True)

    dyn_variables = ["pet", "pre", "tavg"]

    # # loop over dyn variables
    # for var in dyn_variables:
    #     # path for saving aggregated variables
    #     var_new_path = lstm_data_dir + f"dynamic/{var}"
    #     os.makedirs(var_new_path, exist_ok=True)
    #
    #     # path at which netcdf is stored
    #     var_path =  f"/forcings/{var}.nc"
    #
    #     # loop over basins
    #     count = 0
    #     for basin in split_basins:
    #         count += 1
    #         print(f"{split} - {var}: {count}/{len(split_basins)}")
    #         # load netcdf
    #         ds = xr.open_dataset(split_dir + basin + var_path)
    #
    #         # compute mean and save as csv
    #         agg_coord = [coord for coord in list(ds.dims.keys()) if coord != "time"]
    #         mean_var = ds[var].mean(dim=agg_coord)
    #         var_df = mean_var.to_dataframe(name=var).reset_index()
    #         var_df.to_csv(var_new_path + f"/{basin}.csv", sep=";", index=False)
    #         ds.close()

    # Available static variables
    static_variables = ["aspect", "bd", "clay", "dem", "karstic", "lai", "map", "mat", "mat_range", "sand", "slope"]
    # loop over dyn variables
    all_vars = {}
    all_vars["basin"] = split_basins
    for var in static_variables:
        # path at which netcdf is stored
        if var != "dem":
            var_path =  f"/static/mpr/{var}.nc"
        else:
            var_path =  f"/static/routing/{var}.nc"

        var_list = []
        if var == "slope":
            area_list = []
        # loop over basins
        count = 0
        for basin in split_basins:
            count += 1
            print(f"{var}: {count}/{len(split_basins)}")
            # load netcdf
            ds = xr.open_dataset(split_dir + basin + var_path)

            # compute areas from slope netcdf
            if var == "slope":
                dx = ds.lon.diff("lon").mean().item()  # spacing in the x-direction (meters)
                dy = ds.lat.diff("lat").mean().item()  # spacing in the y-direction (meters)
                cell_area_m2 = dx * dy
                cell_area_km2 = cell_area_m2 / 1e6
                valid_cells = ds["slope"].notnull().sum().item()
                total_valid_area_km2 = valid_cells * cell_area_km2
                area_list.append(total_valid_area_km2)
            # compute mean and save as csv
            agg_coord = [dim for dim in ds.dims if dim in ds.coords]
            variable_name = list(ds.data_vars.keys())[0]
            if var != variable_name:
                print(f"{var}: {variable_name}")
            mean_var = ds[variable_name].mean(dim=agg_coord)
            var_list.append(float(mean_var))
            ds.close()

        all_vars[var] = var_list
        if var == "slope":
            all_vars["area"] = area_list



    ##################################
    # Complex derived static variables
    ##################################
    # loop over basins
    count = 0
    aridity = []
    frac_snow = []
    high_prec_freq=[]
    low_prec_freq=[]
    high_prec_dur=[]
    low_prec_dur=[]
    p_season=[]
    frac_forest=[]
    frac_impervious=[]
    frac_pervious=[]
    lai_max = []
    lai_diff = []
    for basin in split_basins:
        # dyn data
        prec = xr.open_dataset(split_dir + basin + "/forcings/pre.nc")
        pet = xr.open_dataset(split_dir + basin + "/forcings/pet.nc")
        tavg = xr.open_dataset(split_dir + basin + "/forcings/tavg.nc")

        # Aridity = mean pet / mean prec; data: prec, pet
        # Compute the mean areal (spatial) values for each day
        agg_coord_prec = [coord for coord in list(prec.dims.keys()) if coord != "time"]
        agg_coord_pet = [coord for coord in list(pet.dims.keys()) if coord != "time"]

        daily_mean_prec = prec.pre.mean(dim=agg_coord_prec)
        daily_mean_pet = pet.pet.mean(dim=agg_coord_pet)
        # Resample to annual totals by summing the daily means within each year
        annual_total_pre = daily_mean_prec.resample(time="1Y").sum(dim="time")
        annual_total_pet = daily_mean_pet.resample(time="1Y").sum(dim="time")
        # Compute the mean annual total by averaging over the years
        mean_annual_pre = annual_total_pre.mean(dim="time")
        mean_annual_pet = annual_total_pet.mean(dim="time")
        basin_aridity = mean_annual_pre / mean_annual_pet
        aridity.append(float(basin_aridity))

        # frac_snow = fraction of prec during days with <0°C, prec, tavg
        snow_days_per_grid = (((prec.pre > 0) & (tavg.tavg < 0)).astype(int)).sum(dim="time")
        # Sum precipitation over time and space for days with tavg < 0°C:
        snow_prec_total = prec.pre.where(tavg.tavg < 0).sum(dim=["time"] + agg_coord_prec)
        # Sum the total precipitation over the entire period and catchment:
        total_prec_total = prec.pre.sum(dim=["time"] + agg_coord_prec)
        # Compute the snow fraction:
        basin_frac_snow = snow_prec_total / total_prec_total
        frac_snow.append(float(basin_frac_snow))

        # high_precip_freq
        mean_daily_prec = daily_mean_prec.mean(dim="time")
        threshold = 5 * mean_daily_prec
        high_precip_days = daily_mean_prec >= threshold
        high_precip_days_per_year = high_precip_days.resample(time="1Y").sum(dim="time")
        basin_high_precip_freq = high_precip_days_per_year.mean(dim="time")
        high_prec_freq.append(float(basin_high_precip_freq))

        # low_precip_freq
        dry_days = daily_mean_prec < 1
        dry_days_per_year = dry_days.resample(time="1Y").sum(dim="time")
        basin_low_precip_freq = dry_days_per_year.mean(dim="time")
        low_prec_freq.append(float(basin_low_precip_freq))

        # high_prec_dur
        # 2. Compute the overall mean daily precipitation.
        mean_daily_prec = daily_mean_prec.mean(dim="time")
        # 3. Define the threshold as 5 times the overall mean daily precipitation.
        threshold = 5 * mean_daily_prec
        # 4. Create a boolean time series indicating high precipitation days.
        high_prec_days = daily_mean_prec >= threshold
        # Convert the boolean DataArray to a numpy array.
        high_prec_np = high_prec_days.values.astype(bool)
        # 5. Define a function to compute the lengths of consecutive True values.
        def compute_run_lengths(arr):
            run_lengths = []
            current_run = 0
            for val in arr:
                if val:
                    current_run += 1
                else:
                    if current_run:
                        run_lengths.append(current_run)
                        current_run = 0
            # Append the last run if the array ends with True.
            if current_run:
                run_lengths.append(current_run)
            return run_lengths
        # 6. Compute the run lengths (i.e., the duration in days for each event).
        event_lengths = compute_run_lengths(high_prec_np)
        # 7. Compute the average duration of high precipitation events.
        if event_lengths:
            avg_duration = np.mean(event_lengths)
        else:
            avg_duration = 0
        high_prec_dur.append(float(avg_duration))

        # low_prec_dur
        # 2. Create a boolean time series indicating dry days (< 1 mm/day).
        dry_days = daily_mean_prec < 1
        # Convert the boolean DataArray to a numpy array.
        dry_days_np = dry_days.values.astype(bool)
        # 3. Function to compute the lengths of consecutive True values.
        def compute_run_lengths(arr):
            run_lengths = []
            current_run = 0
            for val in arr:
                if val:
                    current_run += 1
                else:
                    if current_run:
                        run_lengths.append(current_run)
                        current_run = 0
            # Append the last run if the array ends with True.
            if current_run:
                run_lengths.append(current_run)
            return run_lengths
        # 4. Compute the run lengths (i.e., the duration in days for each dry period).
        dry_period_lengths = compute_run_lengths(dry_days_np)
        # 5. Compute the average duration of dry periods.
        if dry_period_lengths:
            avg_dry_duration = np.mean(dry_period_lengths)
        else:
            avg_dry_duration = 0

        low_prec_dur.append(float(avg_dry_duration))

        # precip_season
        # 1. Compute the daily spatial (areal) mean for precipitation and temperature
        agg_coord_tavg = [coord for coord in list(tavg.dims.keys()) if coord != "time"]
        daily_temp = tavg.tavg.mean(dim=agg_coord_tavg)
        # 2. Extract day-of-year from the time coordinate
        # (xarray has a .dt accessor that provides dayofyear)
        day_of_year = daily_mean_prec.time.dt.dayofyear.values  # numpy array of day-of-year
        # Convert daily time series to numpy arrays for fitting
        prec_values = daily_mean_prec.values
        temp_values = daily_temp.values
        # 3. Compute overall mean values (these serve as the constant term in the sine models)
        P0 = np.nanmean(prec_values)  # mean daily precipitation over the catchment
        T0 = np.nanmean(temp_values)  # mean daily temperature over the catchment
        # 4. Define sine-wave model functions (with fixed mean values)
        def model_temp(t, delta_t, s_t):
            """Temperature model: T0 + delta_t * sin(2*pi*(t - s_t)/365.25)"""
            return T0 + delta_t * np.sin(2 * np.pi * (t - s_t) / 365.25)

        def model_prec(t, delta_p, s_p):
            """Precipitation model: P0 * (1 + delta_p * sin(2*pi*(t - s_p)/365.25))"""
            return P0 * (1 + delta_p * np.sin(2 * np.pi * (t - s_p) / 365.25))
        # 5. Compute an initial guess for s_p (precipitation phase)
        #    Compute mean precipitation for each month and find the month with maximum mean.
        monthly_means = daily_mean_prec.groupby("time.month").mean().values  # shape (12,)
        max_month = np.argmax(monthly_means) + 1  # months are 1-indexed in this context
        s_p_first_guess = 90 - (max_month * 30)
        if s_p_first_guess < 0:
            s_p_first_guess += 360
        # 6. Fit the temperature model to the daily temperature data.
        #    Use initial guesses: delta_t = 5, s_t = -90
        popt_temp, _ = curve_fit(model_temp, day_of_year, temp_values, p0=[5, -90])
        delta_t, s_t = popt_temp
        # 7. Fit the precipitation model to the daily precipitation data.
        #    Use initial guesses: delta_p = 0.4, s_p = s_p_first_guess
        popt_prec, _ = curve_fit(model_prec, day_of_year, prec_values, p0=[0.4, s_p_first_guess])
        delta_p, s_p = popt_prec
        # 8. Compute the p_seasonality variable as defined:
        p_seasonality = delta_p * np.sign(delta_t) * np.cos(2 * np.pi * (s_p - s_t) / 365.25)
        p_season.append(float(p_seasonality))

        # close dyn xarrays
        prec.close()
        pet.close()
        tavg.close()

        # frac fores
        lc = xr.open_dataset(split_dir + basin + "/static/mpr/land_cover.nc")
        # Select the land_cover variable for the last land cover period.
        last_period = lc.land_cover.isel(land_cover_period=-1)

        # Total number of pixels in the last period
        total_pixels = last_period.size

        # Compute the fraction for each class:
        basin_frac_forest = (last_period == 1).sum() / total_pixels
        basin_frac_Impervious = (last_period == 2).sum() / total_pixels
        basin_frac_pervious = (last_period == 3).sum() / total_pixels
        frac_forest.append(float(basin_frac_forest))
        frac_impervious.append(float(basin_frac_Impervious))
        frac_pervious.append(float(basin_frac_pervious))
        lc.close()

        # lai based attributes
        lai = xr.open_dataset(split_dir + basin + "/static/mpr/lai.nc")
        lai_spatial = lai.lai.mean(dim=["lat", "lon"])
        # 2. Compute the maximum monthly mean LAI over the catchment.
        lai_max_catchment = lai_spatial.max(dim="month_of_year")
        # 3. Compute the minimum monthly mean LAI over the catchment.
        lai_min_catchment = lai_spatial.min(dim="month_of_year")
        # 4. Compute lai_diff as the difference between the maximum and minimum monthly means.
        lai_diff_catchment = lai_max_catchment - lai_min_catchment
        lai_max.append(float(lai_max_catchment))
        lai_diff.append(float(lai_diff_catchment))


    # save static data
    all_vars["aridity"] = aridity
    all_vars["frac_snow"] = frac_snow
    all_vars["high_prec_freq"] = high_prec_freq
    all_vars["low_prec_freq"] = low_prec_freq
    all_vars["high_prec_dur"] = high_prec_dur
    all_vars["low_prec_dur"] = low_prec_dur
    all_vars["p_season"] = p_season
    all_vars["frac_forest"] = frac_forest
    all_vars["frac_impervious"] = frac_impervious
    all_vars["frac_pervious"] = frac_pervious
    all_vars["lai_max"] = lai_max
    all_vars["lai_diff"] = lai_diff

    # save static data
    static_data_df = pd.DataFrame(all_vars)
    static_data_df.to_csv(lstm_data_dir + f"static/basin_attributes_{split}.csv", sep=";", index=False)

static_data_1 = pd.read_csv(lstm_data_dir + f"static/basin_attributes_Training.csv", sep=";")
static_data_2 = pd.read_csv(lstm_data_dir + f"static/basin_attributes_Validation.csv", sep=";")
static_data_df = pd.concat([static_data_1, static_data_2], axis=0)
static_data_df = static_data_df.reset_index(drop=True)
static_data_df.to_csv(lstm_data_dir + f"static/basin_attributes.csv", sep=";", index=False)

# runoff data
runoff_dir = path + "/03_discharge/"
basins = [name.split(".")[0] for name in os.listdir(runoff_dir) if not os.path.isdir(os.path.join(runoff_dir, name)) and ".txt" in name]
# path for saving aggregated variables
lstm_data_dir = path + f"/09_lstm_data/"
runoff_new_path = lstm_data_dir + f"dynamic/discharge_m3s"
os.makedirs(runoff_new_path, exist_ok=True)

for basin in basins:
    print(f"{basins.index(basin) + 1}/{len(basins)}")
    header_skip = 0
    with open(runoff_dir + f"{basin}.txt", 'r') as file:
        for i, line in enumerate(file):
            if line.strip().startswith("end"):
                header_skip = i + 1  # skip all rows up to and including the "end" line
                break
    df = pd.read_csv(
        runoff_dir + f"{basin}.txt",
        delim_whitespace=True,  # use any whitespace as delimiter
        skiprows=header_skip,  # skip the first 5 lines (metadata)
        header=None,  # there is no header in the data
        names=['year', 'month', 'day', 'hour', 'minute', 'discharge'],
        na_values=-999  # treat -999 as missing data
    )

    # Combine the date/time columns into a single datetime column.
    df['date'] = pd.to_datetime(df[['year', 'month', 'day', 'hour', 'minute']])
    df = df.drop(columns=['year', 'month', 'day', 'hour', 'minute'])
    df = df.loc[:, ["date", "discharge"]]
    df.to_csv(runoff_new_path + f"/{basin}.csv", sep=";", index=False)


