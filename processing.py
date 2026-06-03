from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyulog import ULog


CP_VALID = 0x02
HALF_CYCLE_VALID = 0x04

GNSS_PREFIX = {
    0: "G",  # GPS
    1: "S",  # SBAS
    2: "E",  # Galileo
    3: "C",  # BeiDou
    5: "J",  # QZSS
    6: "R",  # GLONASS
    7: "I",  # NavIC / IRNSS
}


def valid_flags(flags):
    flags = pd.to_numeric(flags, errors="coerce")
    return (
        flags.notna()
        & ((flags.astype("Int64") & CP_VALID) != 0)
        & ((flags.astype("Int64") & HALF_CYCLE_VALID) != 0)
    )


def rolling_mean_detrend(df, columns, window):
    signal = df[columns].apply(pd.to_numeric, errors="coerce")
    signal = signal.mask(signal == 0)
    rolling_mean = signal.rolling(window=window, center=True, min_periods=1).mean()
    detrended = signal - rolling_mean
    return signal, rolling_mean, detrended


def plot_detrended_on_axis(ax, time, detrended, labels=None):
    for col in detrended.columns:
        label = labels.get(col, col) if labels else col
        ax.plot(time, detrended[col], linewidth=1, label=label)

    ax.set_title("Detrended carrier phase")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Carrier phase - rolling mean (cycles)")
    ax.grid(True, alpha=0.3)

    if len(detrended.columns) <= 16:
        ax.legend(loc="best", fontsize="small")


def get_dataset(ulog, name):
    for dataset in ulog.data_list:
        if dataset.name == name:
            return dataset
    raise ValueError(f"No '{name}' dataset found in the ULG file.")


def satellite_label(gnss_id, sv_id):
    if not np.isfinite(gnss_id) or not np.isfinite(sv_id):
        return ""

    gnss_id = int(gnss_id)
    sv_id = int(sv_id)
    prefix = GNSS_PREFIX.get(gnss_id, f"N{gnss_id}_")
    if len(prefix) == 1:
        return f"{prefix}{sv_id:02d}"
    return f"{prefix}{sv_id}"


def time_seconds(data, t0_us):
    timestamp = data.get("timestamp_sample", data["timestamp"])
    return (np.asarray(timestamp, dtype=np.float64) - t0_us) * 1e-6


def find_time_origin_us(ulog):
    timestamps = []
    for dataset_name in ["vehicle_attitude", "sensor_gps_raw"]:
        data = get_dataset(ulog, dataset_name).data
        timestamp = data.get("timestamp_sample", data["timestamp"])
        timestamps.append(float(np.asarray(timestamp, dtype=np.float64)[0]))
    return min(timestamps)


def quaternion_to_euler_deg(q):
    q0 = q[:, 0]
    q1 = q[:, 1]
    q2 = q[:, 2]
    q3 = q[:, 3]

    roll = np.arctan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2))
    sin_pitch = 2 * (q0 * q2 - q3 * q1)
    pitch = np.arcsin(np.clip(sin_pitch, -1.0, 1.0))
    yaw = np.arctan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3))

    return (
        np.rad2deg(np.unwrap(roll)),
        np.rad2deg(np.unwrap(pitch)),
        np.rad2deg(np.unwrap(yaw)),
    )


def extract_orientation(ulog, t0_us):
    data = get_dataset(ulog, "vehicle_attitude").data
    q_fields = ["q[0]", "q[1]", "q[2]", "q[3]"]
    q = np.column_stack([np.asarray(data[field], dtype=np.float64) for field in q_fields])
    roll_deg, pitch_deg, yaw_deg = quaternion_to_euler_deg(q)

    return pd.DataFrame(
        {
            "time_s": time_seconds(data, t0_us),
            "roll_deg": roll_deg,
            "pitch_deg": pitch_deg,
            "yaw_deg": yaw_deg,
            "q0": q[:, 0],
            "q1": q[:, 1],
            "q2": q[:, 2],
            "q3": q[:, 3],
        }
    )


def optional_field(data, name, length, fill_value=np.nan):
    if name not in data:
        return np.full(length, fill_value)
    return np.asarray(data[name])


def extract_carrier_phase_wide(ulog, t0_us, gps_only=False):
    data = get_dataset(ulog, "sensor_gps_raw").data
    time_s = time_seconds(data, t0_us)
    row_count = len(time_s)

    carrier_fields = sorted(
        field for field in data
        if field.startswith("carrier_phase[") and field.endswith("]")
    )

    rows = []
    for carrier_field in carrier_fields:
        channel = int(carrier_field.split("[", 1)[1].split("]", 1)[0])
        carrier = np.asarray(data[carrier_field], dtype=np.float64)
        gnss_id = optional_field(data, f"gnss_id[{channel}]", row_count)
        sv_id = optional_field(data, f"sv_id[{channel}]", row_count)
        flags = optional_field(data, f"flags[{channel}]", row_count)

        valid_sat = np.isfinite(sv_id) & (sv_id > 0)
        if gps_only:
            valid_sat &= gnss_id == 0
        valid_phase = np.isfinite(carrier) & (carrier != 0)
        valid = valid_sat & valid_phase

        if np.any(np.isfinite(flags)):
            valid &= valid_flags(pd.Series(flags)).to_numpy()

        for index in np.flatnonzero(valid):
            rows.append(
                {
                    "time": time_s[index],
                    "satellite": satellite_label(gnss_id[index], sv_id[index]),
                    "carrier_phase": carrier[index],
                }
            )

    if not rows:
        raise ValueError("No valid carrier-phase samples found in sensor_gps_raw.")

    long = pd.DataFrame(rows)
    wide = (
        long.pivot_table(
            index="time",
            columns="satellite",
            values="carrier_phase",
            aggfunc="mean",
        )
        .sort_index()
        .reset_index()
    )
    wide.columns.name = None
    labels = {col: col for col in wide.columns if col != "time"}
    return wide, labels


def plot_orientation_and_detrended_carrier(orientation, time, detrended, labels, ulg_file):
    fig, (ax_orientation, ax_carrier) = plt.subplots(1, 2, figsize=(16, 7), sharex=False)

    ax_orientation.plot(orientation["time_s"], orientation["roll_deg"], label="roll", linewidth=1.2)
    ax_orientation.set_title("Drone roll")
    ax_orientation.set_xlabel("Time since log start (s)")
    ax_orientation.set_ylabel("Angle (deg)")
    ax_orientation.grid(True, alpha=0.3)
    ax_orientation.legend(loc="best")

    plot_detrended_on_axis(ax_carrier, time, detrended, labels=labels)
    ax_carrier.set_xlabel("Time since log start (s)")

    fig.suptitle(Path(ulg_file).name)
    fig.tight_layout()
    return fig


def run(args):
    ulg_path = Path(args.input_file)
    if ulg_path.suffix.lower() != ".ulg":
        raise ValueError(f"Expected a .ulg file, got: {ulg_path}")

    print(f"Loading ULG file: {ulg_path}")
    ulog = ULog(str(ulg_path))
    t0_us = find_time_origin_us(ulog)

    orientation = extract_orientation(ulog, t0_us)
    carrier_wide, labels = extract_carrier_phase_wide(ulog, t0_us, gps_only=args.gps_only)
    carrier_columns = [col for col in carrier_wide.columns if col != "time"]
    if args.max_satellites is not None:
        carrier_columns = carrier_columns[:args.max_satellites]

    _, rolling_mean, detrended = rolling_mean_detrend(carrier_wide, carrier_columns, args.window)

    print(f"Loaded {len(orientation)} orientation samples from vehicle_attitude.")
    print(f"Detrended {len(carrier_columns)} carrier-phase satellites from sensor_gps_raw.")

    if args.save_prefix:
        prefix = Path(args.save_prefix)
        orientation.to_csv(prefix.with_name(prefix.name + "_orientation.csv"), index=False)
        carrier_wide.to_csv(prefix.with_name(prefix.name + "_carrier_phase_raw_wide.csv"), index=False)
        rolling_mean.to_csv(prefix.with_name(prefix.name + "_carrier_phase_rolling_mean.csv"), index=False)
        detrended.to_csv(prefix.with_name(prefix.name + "_carrier_phase_detrended.csv"), index=False)

    fig = plot_orientation_and_detrended_carrier(
        orientation,
        carrier_wide["time"],
        detrended,
        labels,
        ulg_path,
    )

    if args.save_plot:
        fig.savefig(args.save_plot, dpi=200)
        print(f"Saved plot: {args.save_plot}")

    if not args.no_show:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot drone roll next to detrended GNSS carrier phase from a PX4 ULG file."
    )
    parser.add_argument("input_file", help="Path to input .ulg file.")
    parser.add_argument("--window", type=int, default=50, help="Rolling mean window size in samples.")
    parser.add_argument("--gps-only", action="store_true", help="Only plot GPS satellites.")
    parser.add_argument("--max-satellites", type=int, help="Limit plotted carrier satellites.")
    parser.add_argument("--save-plot", help="Save the combined plot to this image path.")
    parser.add_argument("--save-prefix", help="Save extracted/intermediate CSVs with this prefix.")
    parser.add_argument("--no-show", action="store_true", help="Do not open the plot window.")
    args = parser.parse_args()

    run(args)


if __name__ == "__main__":
    main()
