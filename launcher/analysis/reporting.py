"""CSV reporting helpers for post-simulation analyses."""

import csv
import os


def _format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.10g}"
    return value


def write_csv(path, headers, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({header: _format_value(row.get(header)) for header in headers})
    return path
