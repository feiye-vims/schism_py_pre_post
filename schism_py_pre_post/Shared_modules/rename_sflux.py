"""
Rename sflux file names
Run this script under the run dir or where the old sflux dir is located
"""
import os
import re
from pathlib import Path

# rename the old sflux dir
os.rename('sflux', 'sflux_old')
old_sflux_dir = Path('sflux_old')

os.makedirs('sflux', exist_ok=True)
new_sflux_dir = Path('sflux')

os.chdir(new_sflux_dir)
sflux_nc_files = sorted(Path(f"../{old_sflux_dir}").glob("sflux_*.nc"))

for old_name in sflux_nc_files:
    # Construct new name by stripping leading zeros after first dot
    new_name = re.sub(r'\.(\d+)\.', lambda m: f".{int(m.group(1))}.", old_name.name)

    # Create symbolic link
    os.symlink(f"../{old_sflux_dir / old_name.name}", new_name)
    print(f"Symlink created: {new_name} -> {old_name}")

# link sflux_inputs.txt
os.symlink(f"../{old_sflux_dir / 'sflux_inputs.txt'}", './sflux_inputs.txt')
print(f"Symlink created: {new_sflux_dir / 'sflux_inputs.txt'} -> {old_sflux_dir / 'sflux_inputs.txt'}")