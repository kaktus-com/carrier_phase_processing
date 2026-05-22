# Carrier Phase Processing

This script keeps the same row layout as `flight.csv` while adding minimal carrier phase signals per channel. It converts cycles to meters, computes per-satellite delta range, and skips samples after locktime resets.

## Output fields

The generated `processed_for_plotjuggler.csv` includes all original columns plus new per-channel metrics:

- `phase_m.XX`
- `delta_range_m.XX`

Where `XX` is the channel index (00–31).

## Run

```bash
/home/kate/Documents/Thesis/carrier_phase_processing/.venv/bin/python processing.py --input flight.csv
```

## Options

- `--output`: output CSV path
- `--keep-zero`: keep zero carrier phase samples
- `--time-col`: override time column used for logging (defaults to `sensor_gps_raw/timestamp_sample` if present)

## Notes

- Channels with invalid satellite IDs are ignored when computing delta range.
- Satellites are tracked by (gnss_id, sv_id); if locktime resets, the next measurement is ignored until a new reference is set.
