import base64

from upgrade_modisco_report import available_original_images, resolve_image_bytes


def test_available_original_images_skips_added_columns_and_missing_files(tmp_path):
    payload = b"embedded-image"
    uri = "data:image/png;base64," + base64.b64encode(payload).decode()
    row = {
        "image_records": [
            {"src": str(tmp_path / "missing.png"), "header": "CWM forward"},
            {"src": uri, "header": "PFM forward"},
            {"src": uri, "header": "match0"},
        ]
    }

    images, warnings = available_original_images(row, tmp_path / "report.html", 4)

    assert images == [(payload, "match0")]
    assert len(warnings) == 1
    assert "missing.png" in warnings[0]


def test_available_original_images_applies_limit_after_resolution(tmp_path):
    first = tmp_path / "first.png"
    second = tmp_path / "second.png"
    first.write_bytes(b"first")
    second.write_bytes(b"second")
    row = {
        "image_records": [
            {"src": str(tmp_path / "missing.png"), "header": "CWM forward"},
            {"src": str(first), "header": "CWM forward"},
            {"src": str(second), "header": "CWM reverse"},
        ]
    }

    images, warnings = available_original_images(row, tmp_path / "report.html", 1)

    assert images == [(b"first", "CWM forward")]
    assert len(warnings) == 1


def test_resolve_image_rebases_container_profile_report_path(tmp_path):
    report_dir = tmp_path / "project" / "profile_report"
    logo = report_dir / "trimmed_logos" / "pattern.cwm.fwd.png"
    logo.parent.mkdir(parents=True)
    logo.write_bytes(b"logo")
    report = report_dir / "motifs.html"

    image, warning = resolve_image_bytes(
        "/work/GC_contribs/profile_report/trimmed_logos/pattern.cwm.fwd.png",
        report,
    )

    assert image == b"logo"
    assert warning is None
