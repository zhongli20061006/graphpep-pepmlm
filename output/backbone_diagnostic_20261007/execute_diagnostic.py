"""Execute the approved pilot then full batch, without selecting score direction."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import traceback

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import backbone_diagnostic as diagnostic


def utc():
    return datetime.now(timezone.utc).isoformat()


class LogAndConsole:
    def __init__(self, console, log):
        self.console, self.log = console, log

    def write(self, text):
        self.console.write(text)
        self.log.write(text)

    def flush(self):
        self.console.flush()
        self.log.flush()


def main():
    probe = diagnostic.load_json(OUT / "runtime_preflight/report.json")
    assert probe["status"] == "passed" and probe["excluded_from_formal_statistics"]
    record = dict(started_utc=utc(), status="running", phase="prepare", stages=[],
                  pilot_gate_uses_score_direction=False,
                  implementation_hashes=diagnostic.implementation_hashes(),
                  independent_audit_script_sha256_at_start=diagnostic.sha(OUT / "independent_audit.py"))
    def save():
        diagnostic.save_json(OUT / "run_record.json", record)
    def stage(name, function):
        record["phase"] = name
        item = dict(name=name, started_utc=utc(), status="running")
        record["stages"].append(item)
        save()
        print(f"STAGE {name} started {item['started_utc']}", flush=True)
        result = function()
        item.update(status="completed", finished_utc=utc())
        if name == "prepare":
            result = {k: result[k] for k in ("source_count", "structure_job_count", "pair_count", "contrast_count")}
        item["result"] = result
        save()
        print(json.dumps(dict(stage=name, result=result), ensure_ascii=False, allow_nan=False), flush=True)
        return result
    with (OUT / "execution.log").open("a", encoding="utf-8", buffering=1) as log:
        console = sys.stdout
        sys.stdout = LogAndConsole(console, log)
        try:
            stage("prepare", lambda: diagnostic.prepare(OUT, ROOT))
            stage("baseline", lambda: diagnostic.analysis.baseline(OUT, ROOT, diagnostic.load_json(OUT / "protocol.json")))
            pilot = stage("pilot", lambda: diagnostic.run(OUT, ROOT, "pilot", 2))
            assert pilot["technical_pass"] and pilot["expand_to_full"]
            full = stage("full", lambda: diagnostic.run(OUT, ROOT, "full", 2))
            assert full["technical_pass"] and full["complete"]
            verified = stage("verify", lambda: diagnostic.verify(OUT, ROOT))
            assert verified["complete"]
            record.update(status="computation_verified", finished_utc=utc(), phase="awaiting_independent_audit_and_interpretation")
            save()
        except Exception as error:
            record.update(status="failed", error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
            save()
            raise
        finally:
            sys.stdout = console


if __name__ == "__main__":
    main()
