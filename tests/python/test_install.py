"""AnnotSV-install (port of INSTALL_annotations.sh and the Makefile data targets),
against a local HTTP server: nothing is downloaded from the internet."""
import functools
import http.server
import io
import json
import os
import tarfile
import threading
import zipfile

import pytest

from annotsv import install


def _tar(path, files):
    with tarfile.open(path, "w:gz") as t:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))


class _RangeHandler(http.server.SimpleHTTPRequestHandler):
    """SimpleHTTPRequestHandler + "Range: bytes=N-" (what curl -C - sends)."""
    requests = []

    def log_message(self, *a):
        pass

    def do_GET(self):
        type(self).requests.append((self.path, self.headers.get("Range")))
        rng = self.headers.get("Range")
        path = self.translate_path(self.path)
        if not rng or not os.path.isfile(path):
            return super().do_GET()
        data = open(path, "rb").read()
        start = int(rng.split("=")[1].rstrip("-"))
        if start >= len(data):
            self.send_response(416)
            self.end_headers()
            return
        self.send_response(206)
        self.send_header("Content-Length", str(len(data) - start))
        self.end_headers()
        self.wfile.write(data[start:])


@pytest.fixture
def server(tmp_path, monkeypatch):
    root = tmp_path / "www"
    (root / "ann").mkdir(parents=True)
    (root / "exo").mkdir()
    (root / "jar" / "14.1.0").mkdir(parents=True)
    _tar(root / "ann" / "Annotations_Human_3.5.tar.gz",
         {"Annotations_Human/Genes/GRCh38/genes.txt": b"human\n"})
    _tar(root / "ann" / "Annotations_Mouse_3.4.2.tar.gz",
         {"Annotations_Mouse/Genes/mm39/genes.txt": b"mouse\n"})
    with zipfile.ZipFile(root / "exo" / "2406_phenotype.zip", "w") as z:
        z.writestr("2406_phenotype/2406_phenotype.h2.db", "exomiser")
    (root / "jar" / "14.1.0" / "exomiser-rest-prioritiser-14.1.0.jar").write_bytes(b"0123456789" * 100)

    _RangeHandler.requests = []
    handler = functools.partial(_RangeHandler, directory=str(root))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    monkeypatch.setattr(install, "ANNOTATIONS_URL", f"{base}/ann")
    monkeypatch.setattr(install, "EXOMISER_DATA_URL", f"{base}/exo")
    monkeypatch.setattr(install, "EXOMISER_JAR_URL", f"{base}/jar")
    monkeypatch.setenv("NO_PROXY", "127.0.0.1")
    monkeypatch.setenv("no_proxy", "127.0.0.1")
    yield root
    httpd.shutdown()


def test_human_installs_annotations_and_exomiser(server, tmp_path):
    ann = tmp_path / "AnnotSV_annotations"
    assert install.main(["--annotations-dir", str(ann), "--install-dir", str(tmp_path), "human"]) == 0
    assert (ann / "Annotations_Human" / "Genes" / "GRCh38" / "genes.txt").read_text() == "human\n"
    assert (ann / "Annotations_Exomiser" / "2406" / "2406_phenotype" / "2406_phenotype.h2.db").exists()
    assert not list(ann.glob("*.tar.gz")) and not list(ann.glob("*.zip"))      # archives removed
    assert oct((ann / "Annotations_Human").stat().st_mode & 0o777) == "0o777"


def test_default_annotations_dir_is_install_share(server, tmp_path):
    inst = tmp_path / "prefix"
    assert install.main(["--install-dir", str(inst), "mouse", "exomiser-jar"]) == 0
    share = inst / "share" / "AnnotSV"
    assert (share / "Annotations_Mouse" / "Genes" / "mm39" / "genes.txt").exists()
    jar = share / "jar" / "exomiser-rest-prioritiser-14.1.0.jar"
    assert jar.read_bytes() == b"0123456789" * 100
    # Installed twice: the jar is not downloaded again
    n = len(_RangeHandler.requests)
    assert install.main(["--install-dir", str(inst), "exomiser-jar"]) == 0
    assert len(_RangeHandler.requests) == n


def test_download_resumes_partial_file(server, tmp_path):
    dest = tmp_path / "x.jar"
    dest.write_bytes(b"0123456789" * 30)
    install.download(f"{install.EXOMISER_JAR_URL}/14.1.0/exomiser-rest-prioritiser-14.1.0.jar", str(dest))
    assert dest.read_bytes() == b"0123456789" * 100
    assert _RangeHandler.requests[-1][1] == "bytes=300-"


def test_existing_complete_archive_is_not_downloaded(server, tmp_path):
    ann = tmp_path / "ann"
    ann.mkdir()
    _tar(ann / "Annotations_Mouse_3.4.2.tar.gz", {"Annotations_Mouse/local.txt": b"x"})
    assert install.main(["--annotations-dir", str(ann), "mouse"]) == 0
    assert (ann / "Annotations_Mouse" / "local.txt").exists()
    assert _RangeHandler.requests == []


def test_exomiser_jar_symlink(tmp_path):
    jar = tmp_path / "my.jar"
    jar.write_text("j")
    link = install.install_exomiser_jar(str(tmp_path / "ann"), jar_file=str(jar))
    assert os.path.islink(link) and os.path.basename(link) == "my.jar"


def test_unsafe_tar_member_is_refused(tmp_path):
    _tar(tmp_path / "bad.tar.gz", {"../evil.txt": b"x"})
    with pytest.raises(Exception):
        install._extract_tar(str(tmp_path / "bad.tar.gz"), str(tmp_path / "out"))
    assert not (tmp_path / "evil.txt").exists()


def test_download_error_returns_1(server, tmp_path, capsys):
    assert install.main(["--annotations-dir", str(tmp_path), "--human-version", "9.9", "human"]) == 1
    assert "ERROR" in capsys.readouterr().err


def test_fix_variantconvert_configs(tmp_path):
    cfg = tmp_path / "configs"
    (cfg / "hs1").mkdir(parents=True)
    (cfg / "hs1" / "annotsv3_from_vcf.json").write_text(json.dumps({"h": ["##contig=<ID=chr1>"]}))
    (cfg / "hs1" / "other.json").write_text('"##contig=<ID=chr1>"')
    install.fix_variantconvert_configs(str(cfg))
    assert not (cfg / "hs1").exists()
    assert "ID=1>" in (cfg / "CHM13" / "annotsv3_from_vcf.json").read_text()
    assert "ID=chr1" in (cfg / "CHM13" / "other.json").read_text()
    for genome in ("GRCh37", "GRCh38", "CHM13"):
        assert (cfg / genome / "annotsv3_from_bed.fullsplit.local.json").exists()
