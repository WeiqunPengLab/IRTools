from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import Iterator


BUNDLED_ANNOTATION_ALIASES = {
    "hg38": "gencode.v49.basic.annotation.level1.irtools.gff3.gz",
    "mm39": "gencode.vM38.basic.annotation.level1.irtools.gff3.gz",
}


@contextmanager
def bundled_annotation_path(species: str) -> Iterator[Path]:
    resource_name = BUNDLED_ANNOTATION_ALIASES.get(species, f"{species}_IR_annotation.gtf.gz")
    resource = files("IRTools").joinpath("data", resource_name)
    with as_file(resource) as annotation_path:
        yield annotation_path


@contextmanager
def annotation_input_path(annotation_input: str) -> Iterator[Path]:
    input_path = Path(annotation_input)
    if input_path.exists():
        yield input_path
        return

    data_dir = files("IRTools").joinpath("data")
    alias_name = BUNDLED_ANNOTATION_ALIASES.get(annotation_input)
    candidate_names = (
        alias_name,
        annotation_input,
        f"{annotation_input}.gz",
        f"{annotation_input}_IR_annotation.gtf.gz",
        f"{annotation_input}_IR_annotation.gff3.gz",
        f"{annotation_input}_IR_annotation.gtf",
        f"{annotation_input}_IR_annotation.gff3",
    )

    for candidate_name in candidate_names:
        if candidate_name is None:
            continue
        candidate = data_dir.joinpath(candidate_name)
        if candidate.is_file():
            with as_file(candidate) as annotation_path:
                yield annotation_path
                return

    available = sorted(
        path.name for path in data_dir.iterdir()
        if path.name.endswith((".gtf", ".gtf.gz", ".gff3", ".gff3.gz"))
    )
    raise FileNotFoundError(
        "Annotation input {!r} is neither an existing file nor a bundled IRTools annotation. "
        "Available bundled annotations: {}.".format(
            annotation_input,
            ", ".join(available) if available else "none",
        )
    )
