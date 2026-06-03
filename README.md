# Carrier Phase Processing

Tools for inspecting GNSS carrier phase from two data sources:

- PX4 ULog files (`.ulg`)
- RINEX observation files (`.obs`, `.rnx`, or RINEX text files saved as `.txt`)

The current focus is comparing drone attitude/roll against detrended carrier phase.

## Files

### `processing.py`

Python tool for PX4 `.ulg` files only.

It reads:

- `vehicle_attitude/q[0..3]` for drone attitude
- `sensor_gps_raw/carrier_phase[...]` for GNSS carrier phase
- `sensor_gps_raw/gnss_id[...]` and `sv_id[...]` for satellite labels
- `sensor_gps_raw/flags[...]` to keep valid carrier-phase samples

It plots:

- left: drone roll
- right: detrended carrier phase per satellite

Detrending is done with a centered rolling mean:

```text
detrended = carrier_phase - rolling_mean(carrier_phase)
```

### `processing.m`

MATLAB tool for PX4 `.ulg` files.

It calls `extract_ulg_carrier_phase.py` to decode the ULog, saves a carrier-phase CSV, and plots carrier phase. It can also add high-pass filtered carrier-phase columns, but remember that a 5 Hz high-pass needs data sampled above 10 Hz.

### `processing_t.m`

MATLAB tool for RINEX observation files.

It reads `.obs`, `.rnx`, or `.txt` files containing RINEX observation data and plots `L1C` carrier phase by satellite.

### `extract_ulg_carrier_phase.py`

Python helper used by `processing.m` to extract carrier phase from `.ulg` files into CSV.

## Python ULG Workflow

Run:

```bash
python3 processing.py log_136_2026-5-22-12-04-26.ulg --gps-only
```

Save the plot:

```bash
python3 processing.py log_136_2026-5-22-12-04-26.ulg --gps-only --save-plot log136_roll_carrier.png
```

Save extracted/intermediate CSV files:

```bash
python3 processing.py log_136_2026-5-22-12-04-26.ulg --gps-only --save-prefix log136
```

This writes:

- `log136_orientation.csv`
- `log136_carrier_phase_raw_wide.csv`
- `log136_carrier_phase_rolling_mean.csv`
- `log136_carrier_phase_detrended.csv`

Useful options:

```bash
--window 50          # rolling mean window in samples
--gps-only           # only plot GPS satellites, labels starting with G
--max-satellites N   # limit number of plotted satellites
--no-show            # save/run without opening a plot window
```

Example with a shorter detrending window:

```bash
python3 processing.py log_136_2026-5-22-12-04-26.ulg --gps-only --window 20
```

## MATLAB ULG Workflow

In MATLAB:

```matlab
cd('/home/kate/Documents/Thesis/carrier_phase_processing')
clear processing
rehash toolboxcache

T = processing("log_136_2026-5-22-12-04-26.ulg");
```

Run without plotting:

```matlab
T = processing("log_136_2026-5-22-12-04-26.ulg", "Plot", false);
```

Use a lower high-pass cutoff:

```matlab
T = processing("log_136_2026-5-22-12-04-26.ulg", "HighPassHz", 0.05);
```

Note: the available `sensor_gps_raw` carrier-phase samples in these logs are about 1 Hz. A 5 Hz high-pass cutoff is not possible on 1 Hz data.

## MATLAB RINEX Workflow

Use `processing_t.m` for RINEX files:

```matlab
cd('/home/kate/Documents/Thesis/carrier_phase_processing')
clear processing_t
rehash toolboxcache

T = processing_t("receiver.obs");
```

For a RINEX file saved as text:

```matlab
T = processing_t("my_observations.txt");
```

`processing_t.m` returns a table with:

- `Time`
- `Satellite`
- `CarrierPhaseCycles`
- `CarrierPhaseMeters`

## Satellite Labels

Satellite labels are built from `gnss_id` and `sv_id`.

Common prefixes:

- `G`: GPS
- `S`: SBAS
- `E`: Galileo
- `C`: BeiDou
- `R`: GLONASS
- `J`: QZSS
- `I`: NavIC / IRNSS

For example:

```text
G16  = GPS satellite 16
S123 = SBAS satellite 123
```

Use `--gps-only` in `processing.py` if you only want GPS satellites.

## Notes

- Carrier phase naturally increases for some satellites and decreases for others because it is related to satellite-receiver range rate and Doppler.
- Drone orientation in ULog is stored as a quaternion in `vehicle_attitude/q[0..3]`.
- `processing.py` currently plots roll only, because yaw is not useful for the current antenna-orientation comparison.
- A true sky plot of satellite azimuth/elevation is not available from these ULG files alone. The ULG has satellite IDs, Doppler, carrier phase, and C/N0, but not satellite azimuth/elevation.
