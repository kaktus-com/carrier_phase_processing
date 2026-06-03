from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import pandas as pd


CP_VALID = 0x02
HALF_CYCLE_VALID = 0x04


def load_csv(csv_file):
    print(f"Loading CSV file: {csv_file}")
    return pd.read_csv(csv_file)


def find_flag_columns(df):
    return [col for col in df.columns if "/flags." in col or col.startswith("flags.")]


def carrier_phase_validity(df, flag_columns):
    flags = df[flag_columns].apply(pd.to_numeric, errors="coerce")
    valid_numbers = flags.notna()
    flags = flags.fillna(0).astype("uint16")
    return (((flags & CP_VALID) != 0) & ((flags & HALF_CYCLE_VALID) != 0) & valid_numbers).astype(int)


def plot_validity(time, validity):
    fig, ax = plt.subplots(figsize=(11, 5))
    for col in validity.columns:
        label = col.rsplit("/", 1)[-1]
        ax.step(time, validity[col], where="post", linewidth=1, alpha=0.8, label=label)

    ax.set_title("Carrier phase flag validity")
    ax.set_xlabel("Time")
    ax.set_ylabel("1 = cpValid and halfCyc are both set")
    ax.set_yticks([0, 1])
    ax.set_ylim(-0.1, 1.1)
    ax.grid(True, alpha=0.3)

    if len(validity.columns) <= 12:
        ax.legend(loc="best", fontsize="small")

    fig.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot carrier-phase validity from u-blox flags.")
    parser.add_argument("csv_file", nargs="?", default="flight3.csv")
    args = parser.parse_args()

    df = load_csv(Path(args.csv_file))
    flag_columns = find_flag_columns(df)
    if not flag_columns:
        raise ValueError("No flag columns like 'sensor_gps_raw/flags.00' were found.")

    time = df["time"] if "time" in df.columns else df.index
    validity = carrier_phase_validity(df, flag_columns)
    plot_validity(time, validity)


if __name__ == "__main__":
    main()
