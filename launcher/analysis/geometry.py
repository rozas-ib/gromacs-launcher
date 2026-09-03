"""Geometry helpers shared by structure-based analyses."""


def parse_gro_box_volume_nm3(gro_path):
    """Return the volume of an orthorhombic or triclinic GRO box in nm^3."""
    with open(gro_path, "r", encoding="utf-8", errors="ignore") as gro_file:
        lines = gro_file.readlines()
    if not lines:
        raise ValueError(f"Could not read GRO file '{gro_path}'")

    try:
        values = [float(token) for token in lines[-1].split()]
    except ValueError as exc:
        raise ValueError(f"Invalid GRO box line in '{gro_path}': {lines[-1].strip()}") from exc

    if len(values) == 3:
        volume = values[0] * values[1] * values[2]
    elif len(values) == 9:
        xx, yy, zz, xy, xz, yx, yz, zx, zy = values
        volume = abs(
            xx * (yy * zz - zy * yz)
            - yx * (xy * zz - zy * xz)
            + zx * (xy * yz - yy * xz)
        )
    else:
        raise ValueError(
            f"Unsupported GRO box format in '{gro_path}': expected 3 or 9 values, got {len(values)}"
        )
    if volume <= 0:
        raise ValueError(f"GRO box volume must be positive in '{gro_path}'")
    return volume
