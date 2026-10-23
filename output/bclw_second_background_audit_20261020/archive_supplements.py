"""Retrieve public primary-source supplements and index original Word tables."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parent
SOURCES = {
    "DeBartolo2014_TableS3_Kd.docx": [
        "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1003693.s005"
    ],
    "DeBartolo2014_TableS6_constructs.docx": [
        "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1003693.s008"
    ],
    "Foight2015_supplement1.docx": [
        "https://ars.els-cdn.com/content/image/1-s2.0-S0022283615002995-mmc1.docx",
        "https://pmc.ncbi.nlm.nih.gov/articles/instance/4520770/bin/NIHMS694099-supplement-1.docx"
    ],
    "Foight2015_supplement2.docx": [
        "https://ars.els-cdn.com/content/image/1-s2.0-S0022283615002995-mmc2.docx",
        "https://pmc.ncbi.nlm.nih.gov/articles/instance/4520770/bin/NIHMS694099-supplement-2.docx"
    ],
}


def archive(item: tuple[str, list[str]]) -> dict:
    """Fetch one binary supplement; reject HTML challenge pages as corrupt data."""
    filename, urls = item
    path = ROOT / "raw" / filename
    attempts = []
    for url in urls:
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "peptide-source-audit/1.0"})
            with urllib.request.urlopen(request, timeout=40) as response:
                data = response.read()
                resolved = response.url
            if not data.startswith(b"PK"):
                raise ValueError(f"Expected OOXML zip, got {data[:30]!r}")
            path.write_bytes(data)
            with zipfile.ZipFile(path) as archive_file:
                root = ET.fromstring(archive_file.read("word/document.xml"))
            namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            paragraphs = ["".join(p.itertext()) for p in root.findall("./w:body/w:p", namespace)]
            tables = [[["".join(cell.itertext()) for cell in row.findall("w:tc", namespace)]
                       for row in table.findall("w:tr", namespace)]
                      for table in root.findall(".//w:tbl", namespace)]
            extracted = {"source_file": filename, "paragraphs": paragraphs, "tables": tables}
            (ROOT / "extracted" / f"{path.stem}.json").write_text(
                json.dumps(extracted, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            attempts.append({"url": url, "resolved_url": resolved, "ok": True})
            return {"file": filename, "status": "archived_ooxml", "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(), "attempts": attempts,
                    "retrieved_utc": datetime.now(timezone.utc).isoformat()}
        except (urllib.error.URLError, TimeoutError, ValueError, zipfile.BadZipFile, KeyError) as exc:
            attempts.append({"url": url, "ok": False, "error": str(exc)})
    return {"file": filename, "status": "retrieval_failed", "attempts": attempts}


def main() -> None:
    """Save all retrieval outcomes, including inaccessible files."""
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(archive, SOURCES.items()))
    (ROOT / "supplement_manifest.json").write_text(
        json.dumps({"supplements": results}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for result in results:
        print(result["file"], result["status"], result.get("bytes", 0))


if __name__ == "__main__":
    main()
