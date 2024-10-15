import pandas as pd
import numpy as np

def load_data(file_path):
    """Load data from a CSV file."""
    try:
        data = pd.read_csv(file_path, skiprows=14)
        RefPower = data["Power (W)"] * 1000     # Convert to mW
        x = np.arange(len(RefPower))            # Index for the x-axis
        return x, RefPower
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None
