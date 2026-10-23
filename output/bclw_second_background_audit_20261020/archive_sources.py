"""Archive public primary-source XML and index evidence; never run a model."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
SOURCES = {
    "PMC152889": "Hinds_2003_Bclw_C_terminal",
    "PMC2781529": "Lee_2009_randomized_BH3",
    "PMC3600422": "DeBartolo_2012_binding_models",
    "PMC3631442": "Dutta_2013_Bfl1_ligands",
    "PMC4072508": "DeBartolo_2014_genome_wide",
    "PMC4520770": "Foight_2015_viral_landscape",
    "PMC5464773": "Jenson_2017_PUMA_epistasis",
}
NS = "{http://www.w3.org/1999/xlink}href"
MANIFEST_PATH = ROOT / "source_manifest.json"
PREVIOUS = ({item["pmcid"]: item for item in json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["sources"]}
            if MANIFEST_PATH.exists() else {})


def article_element(data: bytes) -> ET.Element:
    """Accept JATS from Europe PMC or the NCBI E-utilities wrapper."""
    root = ET.fromstring(data)
    article = root if root.tag == "article" else root.find("./article")
    if article is None:
        raise ValueError(f"Expected JATS article, got {root.tag}")
    return article


def flatten(node: ET.Element | None) -> str:
    """Normalize text without dropping a table's labels or footnotes."""
    if node is None:
        return ""
    return re.sub(r"\s+", " ", "".join(node.itertext())).strip()


def archive(item: tuple[str, str]) -> dict:
    """Fetch one source from official public XML endpoints, preserving errors."""
    pmcid, slug = item
    path = ROOT / "raw" / f"{pmcid}.xml"
    attempts: list[dict] = list(PREVIOUS.get(pmcid, {}).get("attempts", []))
    if not path.exists():
        urls = [
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML",
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=pmc&id={pmcid.removeprefix('PMC')}",
        ]
        for url in urls:
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "peptide-source-audit/1.0"})
                with urllib.request.urlopen(request, timeout=45) as response:
                    data = response.read()
                    final_url = response.url
                article_element(data)
                path.write_bytes(data)
                attempts.append({"url": url, "resolved_url": final_url, "ok": True})
                break
            except (urllib.error.URLError, ET.ParseError, ValueError, TimeoutError) as exc:
                attempts.append({"url": url, "ok": False, "error": str(exc)})
    if not path.exists():
        return {"pmcid": pmcid, "slug": slug, "status": "retrieval_failed", "attempts": attempts}
    data = path.read_bytes()
    article = article_element(data)
    title = flatten(article.find("./front/article-meta/title-group/article-title"))
    ids = {node.attrib.get("pub-id-type", "unknown"): flatten(node)
           for node in article.findall("./front/article-meta/article-id")}
    sections = []
    for node in article.findall("./body//sec"):
        heading = flatten(node.find("title"))
        paragraphs = [flatten(p) for p in node.findall("p")]
        sections.append({"heading": heading, "paragraphs": paragraphs})
    tables = []
    for node in article.findall(".//table-wrap"):
        rows = [[flatten(cell) for cell in row if cell.tag in {"td", "th"}]
                for row in node.findall(".//tr")]
        tables.append({"id": node.attrib.get("id"), "label": flatten(node.find("label")),
                       "caption": flatten(node.find("caption")), "rows": rows,
                       "footnotes": flatten(node.find("table-wrap-foot"))})
    supplements = [{"tag": node.tag, "href": node.attrib.get(NS), "text": flatten(node)}
                  for node in article.iter()
                  if node.tag in {"supplementary-material", "media"} and node.attrib.get(NS)]
    text_lines = [title, json.dumps(ids, ensure_ascii=False)]
    for section in sections:
        text_lines += ["\n## " + section["heading"], *section["paragraphs"]]
    for table in tables:
        text_lines += ["\n## " + table["label"], table["caption"],
                       *[" | ".join(row) for row in table["rows"]], table["footnotes"]]
    text_path = ROOT / "extracted" / f"{slug}.txt"
    text_path.write_text("\n".join(text_lines) + "\n", encoding="utf-8")
    index_path = ROOT / "extracted" / f"{slug}.json"
    index_path.write_text(json.dumps({"title": title, "ids": ids, "sections": sections,
                                    "tables": tables, "supplements": supplements},
                                   ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    has_body = article.find("body") is not None
    return {"pmcid": pmcid, "slug": slug, "title": title, "ids": ids,
            "status": "archived_fulltext_jats" if has_body else "metadata_only_jats",
            "has_fulltext_body": has_body, "path": path.relative_to(ROOT).as_posix(),
            "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
            "retrieved_utc": PREVIOUS.get(pmcid, {}).get("retrieved_utc", datetime.now(timezone.utc).isoformat()),
            "attempts": attempts,
            "tables": len(tables), "supplements": supplements}


def main() -> None:
    """Archive bounded sources and save a machine-readable provenance manifest."""
    (ROOT / "raw").mkdir(exist_ok=True)
    (ROOT / "extracted").mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(archive, SOURCES.items()))
    (ROOT / "source_manifest.json").write_text(
        json.dumps({"scope": "bounded primary-source screen; no models or scores", "sources": results},
                   ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for result in results:
        print(result["pmcid"], result["status"], result.get("bytes", 0), result.get("title", ""))


if __name__ == "__main__":
    main()
