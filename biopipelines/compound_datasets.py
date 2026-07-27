# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Public compound datasets usable as a CompoundLibrary source.

Each dataset is a CompoundDataset descriptor naming a remote archive and the
member CSV to extract. Resolution downloads the archive once into the cache
folder, then normalizes the member CSV to the compounds schema.
"""

import csv
import hashlib
import os
import shutil
import urllib.request
import zipfile
from typing import Dict, List, Optional

# One row per fluorophore-solvent pair in every FluoDB-family table, so the same
# molecule recurs across solvents; dedupe on smiles for a per-molecule library.
_FLUODB_ARCHIVE = "https://ndownloader.figshare.com/files/51101936"
_FLUODB_MD5 = "b7d41872f282ae99077a126cd2b84bfb"
_FLUODB_CITATION = (
    "Wang et al., A modular artificial intelligence framework to facilitate "
    "fluorophore design. Nat Commun 16 (2025). doi:10.1038/s41467-025-58881-5. "
    "Dataset: doi:10.6084/m9.figshare.26317933 (CC BY 4.0)."
)

# Source column -> compounds-schema column. Sub-databases capitalize differently
# and omit varying subsets, so each name is mapped explicitly.
_COLUMN_ALIASES = {
    "smiles": "smiles",
    "absorption/nm": "absorption_nm",
    "emission/nm": "emission_nm",
    "plqy": "plqy",
    "e/m-1cm-1": "extinction",
    "solvent": "solvent",
    "reference(doi)": "reference_doi",
    "source": "source",
    "tag_name": "tag_name",
}

_EXTRA_COLUMNS = [
    "absorption_nm", "emission_nm", "plqy", "extinction",
    "solvent", "reference_doi", "source", "tag_name",
]


class CompoundDataset:
    """A downloadable compound dataset resolvable to a normalized CSV."""

    def __init__(self, name: str, archive_url: str, member: str,
                 archive_md5: Optional[str] = None, citation: str = "",
                 license: str = "", description: str = ""):
        self.name = name
        self.archive_url = archive_url
        self.member = member
        self.archive_md5 = archive_md5
        self.citation = citation
        self.license = license
        self.description = description

    def __repr__(self) -> str:
        return f"CompoundDataset({self.name})"

    def resolve(self, cache_folder: str) -> str:
        """Return the path to this dataset's normalized CSV, downloading once if needed."""
        dataset_dir = os.path.join(cache_folder, "CompoundDatasets")
        normalized = os.path.join(dataset_dir, f"{self.name}.csv")
        if os.path.exists(normalized):
            return normalized

        os.makedirs(dataset_dir, exist_ok=True)
        archive = os.path.join(dataset_dir, f"_{_archive_key(self.archive_url)}.zip")
        if not os.path.exists(archive):
            _download(self.archive_url, archive, self.archive_md5)

        with zipfile.ZipFile(archive) as zf:
            if self.member not in zf.namelist():
                raise ValueError(
                    f"Dataset '{self.name}': member '{self.member}' not found in "
                    f"{self.archive_url}. The upstream archive layout may have changed."
                )
            with zf.open(self.member) as src:
                rows, fieldnames = _read_csv(src)

        _write_normalized(normalized, rows, fieldnames, self.name)
        return normalized


def _archive_key(url: str) -> str:
    return hashlib.md5(url.encode()).hexdigest()[:12]


def _download(url: str, destination: str, expected_md5: Optional[str]):
    """Download to a temporary path, verify, then move into place."""
    partial = destination + ".partial"
    print(f"  Downloading compound dataset archive: {url}")
    with urllib.request.urlopen(url) as response, open(partial, "wb") as out:
        shutil.copyfileobj(response, out)

    if expected_md5:
        digest = hashlib.md5()
        with open(partial, "rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        if digest.hexdigest() != expected_md5:
            os.remove(partial)
            raise ValueError(
                f"Checksum mismatch for {url}: expected {expected_md5}, "
                f"got {digest.hexdigest()}. Refusing to use the download."
            )

    os.replace(partial, destination)


def _read_csv(handle):
    import io
    text = io.TextIOWrapper(handle, encoding="utf-8-sig", errors="replace", newline="")
    reader = csv.DictReader(text)
    rows = list(reader)
    return rows, list(reader.fieldnames or [])


def _write_normalized(path: str, rows: List[Dict], fieldnames: List[str], dataset_name: str):
    """Map source columns onto the compounds schema and write the cached CSV."""
    mapping = {}
    for column in fieldnames:
        if not column or column.startswith("Unnamed:"):
            continue
        target = _COLUMN_ALIASES.get(column.strip().lower())
        if target:
            mapping[column] = target

    if "smiles" not in mapping.values():
        raise ValueError(
            f"Dataset '{dataset_name}': no SMILES column among {fieldnames}. "
            "The upstream schema may have changed."
        )

    header = ["id", "format", "smiles", "ccd"] + _EXTRA_COLUMNS
    written = 0
    with open(path + ".partial", "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(header)
        for index, row in enumerate(rows):
            record = {target: (row.get(source) or "").strip()
                      for source, target in mapping.items()}
            smiles = record.get("smiles", "")
            if not smiles:
                continue
            writer.writerow(
                [f"{dataset_name}_{index}", "smiles", smiles, ""]
                + [record.get(column, "") for column in _EXTRA_COLUMNS]
            )
            written += 1

    if not written:
        os.remove(path + ".partial")
        raise ValueError(f"Dataset '{dataset_name}': no rows carried a SMILES value.")

    os.replace(path + ".partial", path)
    print(f"  Compound dataset '{dataset_name}': {written} compounds cached at {path}")


def _fluodb(name: str, member: str, description: str) -> CompoundDataset:
    return CompoundDataset(
        name=name,
        archive_url=_FLUODB_ARCHIVE,
        member=member,
        archive_md5=_FLUODB_MD5,
        citation=_FLUODB_CITATION,
        license="CC BY 4.0",
        description=description,
    )


FLUODB = _fluodb(
    "FluoDB", "data/FluoDB-Full.csv",
    "Aggregated fluorophore database: 55169 fluorophore-solvent pairs, 35535 unique SMILES.")

FLUODB_LITE = _fluodb(
    "FluoDB-Lite", "data/FluoDB-Lite.csv",
    "FluoDB with complexes and mixtures removed: 49851 pairs, 30939 unique SMILES.")

CHEMFLUOR = _fluodb(
    "ChemFluor", "data/Database/ChemFluor_processed.csv",
    "ChemFluor source database (4052 rows).")

DEEP4CHEM = _fluodb(
    "Deep4Chem", "data/Database/Deep4Chem_processed.csv",
    "Deep4Chem source database (20014 rows).")

DYES = _fluodb(
    "DYES", "data/Database/DYES_processed.csv",
    "DYES source database (25896 rows); absorption and emission absent.")

PHOTOCHEMCAD = _fluodb(
    "PhotochemCAD", "data/Database/PhotochemCAD_processed.csv",
    "PhotochemCAD source database (427 rows); emission absent.")

DSSCDB = _fluodb(
    "DSSCDB", "data/Database/DSSCDB_processed.csv",
    "Dye-sensitized solar cell database (2553 rows).")

CHEMDATAEXTRACTOR = _fluodb(
    "ChemDataExtractor", "data/Database/ChemDataExtractor_processed.csv",
    "ChemDataExtractor-mined source database (2219 rows); emission and PLQY absent.")

DYE_AGGREGATION = _fluodb(
    "DyeAggregation", "data/Database/Dye aggregation_processed.csv",
    "Dye aggregation source database (4043 rows); emission absent.")
