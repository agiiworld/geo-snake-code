# GEOgik <3

A student geodesy project that allows performing basic geodetic calculations and coordinate transformations in Python.

## Features

When you run the program, you’ll see a menu with the following options:

1. **Quasigeoid separation**
   - Bilinear interpolation based on the `WAWA.txt` file.
   - Input: geographic coordinates (latitude and longitude in decimal degrees, with `.` as separator).
   - The program finds the 4 nearest grid points and performs interpolation.
   - Output: separation between quasigeoid and ellipsoid in meters.

2. **Transformation φ, λ → PL-1992**
   - Converts geographic coordinates to rectangular coordinates in the PL-1992 system.

3. **Transformation φ, λ → PL-2000**
   - Converts geographic coordinates to rectangular coordinates in the PL-2000 system.
   - Requires providing the central meridian (`λ0`).

4. **Degrees → Radians**
   - Simple unit conversion.

5. **Radians → Degrees**
   - Result returned in degrees, minutes, and seconds (DMS).

6. **Transformation XYZ → φ, λ, h (Hirvonen’s algorithm)**
   - Converts Cartesian coordinates to ellipsoidal coordinates.

7. **DMS → Decimal degrees**
   - Converts an angle in degrees, minutes, and seconds to decimal degrees.

---

## Input data

- Coordinates should be entered in **decimal degrees**, e.g.:
  ```text
  52.098 21.032
  ```
- If you have data in **degrees, minutes, seconds**, use option **7** to convert them to decimal degrees.

---

## Requirements

- Python 3.x
- Libraries:
  - `numpy`
  - `math`
  - `argparse`

The `WAWA.txt` file (quasigeoid model) **must be in the same folder as the program**.  
➡️ Do not delete or rename it!

---

## How to run

Simply run the program without parameters:

```bash
python proj2_funkcje.py
```

This will open the main menu where you can choose an option.

To display help:

```bash
python proj2_funkcje.py --help
```

---

## Authors

Project created as part of geodesy coursework. 💚
