import pandas as pd
import numpy as np


class DataParser:
    def parse_to_np(path: str) -> np.ndarray:
        df = pd.read_csv(path, sep="\t", skiprows=1)
        df.columns = ["X", "Y"]
        arr = df.to_numpy()
        return arr