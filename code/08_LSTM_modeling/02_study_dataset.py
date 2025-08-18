from pathlib import Path
from typing import List, Dict, Union
import pandas as pd
import xarray
from neuralhydrology.datasetzoo.basedataset import BaseDataset
from neuralhydrology.utils.config import Config
import os
from neuralhydrology.datasetzoo import register_dataset
data_dir = Path("C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/09_lstm_data")

# helper function for normalizing discharge
def _normalize_discharge(ser: pd.Series, area: float) -> pd.Series:
    """Helper function to normalize discharge data by basin area"""
    return ser / (area * 1e6) * 1000 * 86400

# preprocess dynamic data
dyn_dir = Path(data_dir / "dynamic")
preprocessed_dir = Path(data_dir / "preprocessed")
os.makedirs(preprocessed_dir, exist_ok=True)
attributes = pd.read_csv(data_dir / "static" / "basin_attributes.csv", sep=";")

basins = []
for split in ["Training", "Validation"]:
    split_dir = f"C:/Users/morit/Dropbox/Projekte/GenAI_para_for_HLSMs/data/06_training_validation_data/{split}/"
    split_basins = [name for name in os.listdir(split_dir) if os.path.isdir(os.path.join(split_dir, name))]
    basins =  basins + split_basins
count = 0
for basin in basins:
    count += 1
    print(f"{count}/{len(basins)}\n")
    basin_data = []
    for var in ['discharge_m3s', 'pet', 'pre', 'tavg']:
        var_data = pd.read_csv(dyn_dir / var / (basin + ".csv"), sep=";")
        var_data.columns = ["date", var]

        if var == "discharge_m3s":
            basin_area = attributes["area"].loc[[str(basin_id) in basin for basin_id in attributes["basin"]]].item()
            var_data["discharge"] = _normalize_discharge(var_data["discharge_m3s"], basin_area)
            var_data = var_data.drop(columns=["discharge_m3s"])

        var_data['date'] = pd.to_datetime(var_data['date'])
        var_data = var_data.set_index('date')
        basin_data.append(var_data)

    # Merge all variables
    basin_data = pd.concat(basin_data, axis=1)

    # Create complete date range
    full_range = pd.date_range(start=basin_data.index.min(), end=basin_data.index.max(), freq='D')

    # Reindex to include all dates, fill missing with NaN
    basin_data = basin_data.reindex(full_range)
    basin_data.index.name = "date"  # Optional, for clarity in the output CSV

    basin_data.to_csv(preprocessed_dir / (basin + ".csv"))


# from neuralhydrology.datautils.utils import infer_frequency
def load_GenAI_para_ger_timeseries(data_dir: Path, basin: str) -> pd.DataFrame:
    preprocessed_dir = data_dir / "preprocessed"

    # load the data for the specific basin into a time-indexed dataframe
    basin_file = preprocessed_dir / f"{basin}.csv"
    df = pd.read_csv(basin_file, index_col='date', parse_dates=['date'])
    return df



def load_GenAI_para_ger_attributes(data_dir: Path, basins: List[str] = []) -> pd.DataFrame:

    # load attributes into basin-indexed dataframe
    df = pd.read_csv(data_dir / "static" / "basin_attributes.csv", index_col="basin", sep=";")
    df.index = df.index.astype(str)
    # convert all columns, where possible, to numeric
    df = df.apply(pd.to_numeric)

    if basins:
        if any(b not in df.index for b in basins):
            raise ValueError('Some basins are missing static attributes.')
        df = df.loc[basins]

    return df

class GenAI_para_GER(BaseDataset):

    def __init__(self,
                 cfg: Config,
                 is_train: bool,
                 period: str,
                 basin: str = None,
                 additional_features: List[Dict[str, pd.DataFrame]] = [],
                 id_to_int: Dict[str, int] = {},
                 scaler: Dict[str, Union[pd.Series, xarray.DataArray]] = {}):

        # Initialize `BaseDataset` class
        super(GenAI_para_GER, self).__init__(cfg=cfg,
                                       is_train=is_train,
                                       period=period,
                                       basin=basin,
                                       additional_features=additional_features,
                                       id_to_int=id_to_int,
                                       scaler=scaler)

    def _load_basin_data(self, basin: str) -> pd.DataFrame:
        """Load timeseries data of one specific basin"""
        return load_GenAI_para_ger_timeseries(data_dir=self.cfg.data_dir, basin=basin)

    def _load_attributes(self) -> pd.DataFrame:
        """Load catchment attributes"""
        return load_GenAI_para_ger_attributes(self.cfg.data_dir, basins=self.basins)

# register new dataset
register_dataset("GenAI_para_ger", GenAI_para_GER)