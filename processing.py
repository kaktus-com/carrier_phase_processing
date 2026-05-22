import pandas as pd
import numpy as np
from scipy import signal

def parse_carrier_phase(flight):
    print(f"Loading CSV file: {flight}")
    df = pd.read_csv(flight)
    print(f"Total rows: {len(df)}")
    return df

def carrier_phase_processing(col_data, time_data):
    λ = 0.190293672798
    
    # keep structure
    col_data = col_data.copy()
    
    # 1. treat only real measurements as valid
    col_data[col_data == 0] = np.nan
    print(col_data.to_string())
    
    # 2. convert to meters
    phase_meters = col_data * λ
    print(phase_meters.to_string())

    # 2b. high-pass filter (20 Hz cutoff)
    dt = time_data.diff()
    median_dt = dt.dropna().median()
    phase_hp = pd.Series(np.nan, index=phase_meters.index)
    if pd.notna(median_dt) and median_dt > 0:
        fs = 1.0 / median_dt
        cutoff_hz = 1.0
        if fs > 2 * cutoff_hz:
            b, a = signal.butter(4, cutoff_hz / (0.5 * fs), btype="high")
            phase_filled = phase_meters.interpolate(limit_direction="both")
            phase_hp = pd.Series(signal.filtfilt(b, a, phase_filled), index=phase_meters.index)
    
    # 3. diff only where valid
    delta = phase_meters.diff()
    print(delta.to_string())
    dt = time_data.diff()
    
    # 4. rate
    phase_rate = delta
    
    # NORMALIZE RATE: subtract first valid rate to center at zero
    first_valid_rate = phase_rate.dropna().iloc[0] if len(phase_rate.dropna()) > 0 else 0
    phase_rate = phase_rate - first_valid_rate
    print(phase_rate.to_string())
    
    # 5. remove insane spikes (optional but important)
    # phase_rate[phase_rate.abs() > 100] = np.nan
    
    return phase_rate, phase_meters, phase_hp

def main():
    flight = "flight3.csv"
    df_raw = parse_carrier_phase(flight)
    
    if len(df_raw) == 0:
        print("\nNo data loaded from CSV")
        return
    
    # Find all carrier phase columns
    carrier_phase_cols = [
        col for col in df_raw.columns
        if "carrier_phase." in col
    ]
    
    print(f"\nFound {len(carrier_phase_cols)} carrier phase columns")
    
    results = pd.DataFrame()
    
    # Keep time column
    results["time"] = df_raw["time"]
    
    # Process each carrier phase column
    for col in carrier_phase_cols:
        print(f"Processing {col}")
        phase_rate, phase_meters, phase_hp = carrier_phase_processing(
            df_raw[col],
            df_raw["time"]
        )
        results[f"{col}_rate"] = phase_rate
        results[f"{col}_meters"] = phase_meters
        results[f"{col}_meters_hp20"] = phase_hp
    
    # Save
    results.to_csv(
        "processed_carrier_phase3.3.csv",
        index=False
    )
    print("\nSaved to processed_carrier_phase.csv")

if __name__ == "__main__":
    main()