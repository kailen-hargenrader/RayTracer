## Testing the Blender Import/Export Pipeline

This explains how to verify that importing a scene from JSON and exporting it back to JSON produces identical files.

### Prerequisites
- Blender installed and available on your PATH (or know its full path)
- This repository checked out locally

### 1) Import a test scene into Blender
- Use one of the JSON files in `tests/export/` as input.
- In Blender, run the import script `tests/generate_from_json.py` and pass the chosen JSON file.
- After the script runs, visually verify in Blender that the scene looks correct (objects, transforms, materials, cameras, lights as applicable).

Optional command-line example (PowerShell):

```powershell
blender --factory-startup --python tests\generate_from_json.py -- tests\export\<input-scene>.json
```

Notes:
- Arguments after `--` are passed to the Python script.
- Replace `<input-scene>.json` with the actual file name in `tests\export\`.

### 2) Export the scene back to JSON
- From the same Blender session (with the verified scene loaded), run the exporter at `s2899221/Blend/Export` to write out a new JSON file (e.g., `tests/export/<input-scene>.out.json`).

Optional command-line example (PowerShell):

```powershell
blender --factory-startup --python s2899221\Blend\Export.py -- --output tests\export\<input-scene>.out.json
```

Notes:
- If the exporter takes different or additional arguments, use those as required by the script. The general pattern is to pass script-specific args after `--`.

### 3) Verify input and output JSON are identical
If both scripts work correctly, the original input JSON and the exported output JSON should be identical.

Quick comparison options (PowerShell):

```powershell
# No output means files are identical
fc tests\export\<input-scene>.json tests\export\<input-scene>.out.json
```

Or using Python for a structural (whitespace-insensitive) check:

```powershell
python - << 'PY'
import json, sys, pathlib
inp = pathlib.Path(r"tests/export/<input-scene>.json")
out = pathlib.Path(r"tests/export/<input-scene>.out.json")
equal = json.load(open(inp)) == json.load(open(out))
print("IDENTICAL" if equal else "DIFFERENT")
sys.exit(0 if equal else 1)
PY
```

Outcome:
- IDENTICAL/no differences: test passed
- DIFFERENT/differences shown: investigate the import/export scripts or scene contents

