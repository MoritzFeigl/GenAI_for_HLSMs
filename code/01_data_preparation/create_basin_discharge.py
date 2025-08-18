"""
This script generates the discharge data for all basins based on data from the GRDC so that the resulting text file
matches the format of the original mHM discharge input (provided by Moritz).
The raw data can be downloaded from the GRDC download center in the GRDC Export format (ASCII text) and the main
directory containing the text files should be specified as data_dir.
"""
from pathlib import Path

import pandas as pd


def get_attributes(content_list, name_start):
    attr = [line for line in content_list if line.startswith(name_start)][0]
    idx = attr.rfind('  ') + 2
    return attr[idx:-1]


if __name__ == '__main__':

    base_dir = Path('/data/CFG')
    data_dir = base_dir / 'required_data' / 'discharge_GRDC'
    output_dir = base_dir / 'discharge'
    basins_small = pd.read_csv(base_dir / 'required_data' / 'basins' / 'small_basins.csv')
    basins_large = pd.read_csv(base_dir / 'required_data' / 'basins' / 'large_basins.csv')
    basins = pd.concat([basins_large, basins_small])

    files = list(data_dir.glob('*.txt'))
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # List all sites in the GRDC data
    df = basins[['Stat_ID', 'Station_Name', 'River']]
    df = df.set_index('Stat_ID', drop=False)
    df = df[~df.index.duplicated(keep='last')]

    df_grdc = pd.DataFrame(columns=['Station_Name_GRDC', 'River_GRDC'])

    # Process data with matching GRDC ID and Stat_ID
    for n, file in enumerate(files):
        with open(file, encoding='latin1') as f:
            content = f.readlines()

        grdc_id = int(get_attributes(content_list=content, name_start='# GRDC-No.:'))

        river = get_attributes(content_list=content, name_start='# River:')
        station = get_attributes(content_list=content, name_start='# Station:')

        df_grdc.loc[grdc_id, 'Station_Name_GRDC'] = station
        df_grdc.loc[grdc_id, 'River_GRDC'] = river

    df_grdc['Stat_ID_GRDC'] = df_grdc.index

    # good_basins have the same IDs
    good_basins = df_grdc.merge(df, how='inner', left_index=True, right_index=True)

    other_basins = df_grdc.loc[~df_grdc.index.isin(good_basins.index), :].copy()

    # wrong_id_basins have different IDs but share the exact same river and station names
    wrong_id_basins = other_basins.merge(df,
                                         how='inner',
                                         left_on=['Station_Name_GRDC', 'River_GRDC'],
                                         right_on=['Station_Name', 'River'])
    wrong_id_basins = wrong_id_basins[~wrong_id_basins['Stat_ID_GRDC'].duplicated(keep='last')]
    wrong_id_basins = wrong_id_basins.set_index('Stat_ID', drop=False)

    other_basins = other_basins.loc[~other_basins.index.isin(wrong_id_basins['Stat_ID_GRDC']), :]
    other_basins['Station_Name_GRDC'] = other_basins['Station_Name_GRDC'].apply(lambda x: x.replace('-', ' '))
    other_basins['Station_Name_GRDC'] = other_basins['Station_Name_GRDC'].apply(lambda x: x.replace(' 1', ''))
    other_basins['River_GRDC'] = other_basins['River_GRDC'].apply(lambda x: x.replace(' RIVER', ''))

    df_remaining = df.loc[~df.index.isin(good_basins.index) & ~df.index.isin(wrong_id_basins.index), :].copy()
    df_remaining['River'] = df_remaining['River'].apply(lambda x: x.replace('-', ' '))
    df_remaining['Station_Name'] = df_remaining['Station_Name'].apply(lambda x: x.replace('-', ' '))

    # wrong_name_basins have different IDs and some small differences in the names, e.g. '-' and ' '
    wrong_name_basins = other_basins.merge(df_remaining,
                                           how='inner',
                                           left_on=['Station_Name_GRDC', 'River_GRDC'],
                                           right_on=['Station_Name', 'River'])
    wrong_name_basins = wrong_name_basins.set_index('Stat_ID', drop=False)

    df_remaining = df_remaining.loc[~df_remaining.index.isin(wrong_name_basins.index), :].copy()
    other_basins = other_basins.loc[~other_basins.index.isin(wrong_name_basins['Stat_ID_GRDC']), :].copy()

    wrong_basins = other_basins.copy()

    cols = ['Stat_ID', 'Station_Name', 'River']

    # Kirchenhausen
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'KIRCHEN HAUSEN') & (wrong_basins['River_GRDC'] == 'DANUBE')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'KIRCHENHAUSEN') &
            (df_remaining['River'] == 'DONAU'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Gutach
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'GUTACH ELZ') & (wrong_basins['River_GRDC'] == 'ELZ')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'GUTACH') &
            (df_remaining['River'] == 'ELZ'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Geislingen
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'GEISLINGEN FILS') & (wrong_basins['River_GRDC'] == 'FILS')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'GEISLINGEN (BRUECKE)') &
            (df_remaining['River'] == 'FILS'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Grossdittmannsdorf
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'GROSSDITTMANNSDORF') & (
            wrong_basins['River_GRDC'] == 'GROSSE RODER')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'GROSSDITTMANNSDORF') &
            (df_remaining['River'] == 'GROSSE ROEDER'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Goettingen
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'GOETTINGEN') & (wrong_basins['River_GRDC'] == 'LEINE')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'GOTTINGEN') &
            (df_remaining['River'] == 'LEINE'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Dohna
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'DOHNA') & (wrong_basins['River_GRDC'] == 'MUGLITZ')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'DOHNA') &
            (df_remaining['River'] == 'MUEGLITZ'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Rotenfels
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'BAD ROTENFELS') & (wrong_basins['River_GRDC'] == 'MURG')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'ROTENFELS') &
            (df_remaining['River'] == 'MURG'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Streckenwalde
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'STRECKEWALDE') & (wrong_basins['River_GRDC'] == 'PRESSNITZ')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'STRECKENWALDE') &
            (df_remaining['River'] == 'PRESSNITZ'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Rappelsdorf
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'RAPPELSDORF') & (wrong_basins['River_GRDC'] == 'SCHLEUSE')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'RAPPELDORF') &
            (df_remaining['River'] == 'SCHLEUSE'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Stein
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'STEIN') & (wrong_basins['River_GRDC'] == 'TRAUN')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'STEIN BEI ALTENMARKT') &
            (df_remaining['River'] == 'TRAUN'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    # Harburg
    left_condition = (wrong_basins['Station_Name_GRDC'] == 'HARBURG') & (wrong_basins['River_GRDC'] == 'WORNITZ')
    right = df_remaining.loc[
            (df_remaining['Station_Name'] == 'HARBURG') &
            (df_remaining['River'] == 'WOERNITZ'), :]
    wrong_basins.loc[left_condition, cols] = right[cols].values

    wrong_basins = wrong_basins.dropna(axis=0, how='any')
    wrong_basins['Stat_ID'] = wrong_basins['Stat_ID'].astype(int)
    wrong_basins = wrong_basins.set_index('Stat_ID', drop=False)

    available_basins = pd.concat([good_basins, wrong_id_basins, wrong_name_basins, wrong_basins])

    # Adorf was replaced by new gauge during training period (will be discarded)
    adorf_id = int(available_basins.loc[available_basins['Stat_ID'].duplicated(), 'Stat_ID'])
    available_basins = available_basins.loc[available_basins['Stat_ID'] != adorf_id, :]

    # Process data with matching GRDC ID and Stat_ID
    for grdc_id in available_basins['Stat_ID_GRDC'].tolist():

        file = [file for file in files if str(grdc_id) in file.stem][0]
        with open(file, encoding='latin1') as f:
            content = f.readlines()

        river = available_basins.loc[available_basins['Stat_ID_GRDC'] == grdc_id, 'River']
        station = available_basins.loc[available_basins['Stat_ID_GRDC'] == grdc_id, 'Station_Name']
        write_id = int(available_basins.loc[available_basins['Stat_ID_GRDC'] == grdc_id, 'Stat_ID'])

        q_start_index = content.index('# DATA\n') + 2
        dates = []
        values = []
        for line in content[q_start_index:]:
            dates.append(line[:10])
            values_idx = line.rfind('  ') + 2
            values.append(float(line[values_idx:-1]))

        dates = pd.to_datetime(dates)

        available_basins.loc[write_id, 'Start'] = dates.min()
        available_basins.loc[write_id, 'End'] = dates.max()

        first_line = f'     {write_id} {station}-{river}(daily discharge)'
        second_line = f'nodata   -999.'
        third_line = f'n       1       measurements per day [1, 1440]'

        year = dates.min().year
        month = str(dates.min().month).zfill(2)
        day = str(dates.min().day).zfill(2)
        forth_line = f'start  {year} {month} {day} 00 00   (YYYY MM DD HH MM)'

        year = dates.max().year
        month = str(dates.max().month).zfill(2)
        day = str(dates.max().day).zfill(2)
        fifth_line = f'end    {year} {month} {day} 00 00   (YYYY MM DD HH MM)'

        with open(output_dir / f'{write_id}.txt', 'w') as f:
            for line in [first_line, second_line, third_line, forth_line, fifth_line]:
                f.write(line)
                f.write('\n')
        f.close()

        with open(output_dir / f'{write_id}.txt', 'a') as f:
            for (date, value) in zip(dates, values):
                date_str = f'{date.year}  {str(date.month).zfill(2)}  {str(date.day).zfill(2)}  00  00'
                value_str = f'{value:.3f}'
                white_spaces = ' '*(7 - len(str(int(value))))
                f.write(f'{date_str}{white_spaces}{value_str}')
                f.write('\n')
        f.close()

    no_period = available_basins.loc[available_basins['End'] < pd.Timestamp('2014-12-31'), :]
    full_period = available_basins.loc[available_basins['End'] >= pd.Timestamp('2019-12-31'), :]

    available_basins.to_csv(output_dir / 'discharge_gauges_overview.csv', index_label='ID')


