"""Read-only inventory of the existing OpenMM engine and force-field resources."""
import hashlib
import json
from pathlib import Path

import openmm
import openmm._openmm
import openmm.app
import openmm.app.forcefield
import openmm.app.internal.customgbforces

app_dir = Path(openmm.app.__file__).parent
files = [app_dir / "data/amber14/protein.ff14SB.xml", app_dir / "data/implicit/obc2.xml",
         Path(openmm.__file__), Path(openmm._openmm.__file__), Path(openmm.app.forcefield.__file__),
         Path(openmm.app.internal.customgbforces.__file__)]
hashes = {}
for path in files:
    with path.open("rb") as stream:
        hashes[str(path)] = hashlib.file_digest(stream, "sha256").hexdigest()
platform = openmm.Platform.getPlatformByName("CPU")
print(json.dumps({"openmm_version": openmm.__version__, "platform": "CPU", "platform_properties": list(platform.getPropertyNames()),
                  "python_executable": __import__("sys").executable, "external_hashes": hashes}, indent=2))
