from ..typing import NamedTuple, List, Optional


class Sequence(NamedTuple):
    """Sequence which contains the index of the monomer in the monomers list, the information if it is inverted and the number of repetitions."""

    index: List[int]
    inverted: List[bool]
    reps: List[int]


class Atom(NamedTuple):
    """A Namedtuple which contains the information about the atom."""

    ff_identifier: Optional[str]
    element: str
    x: float
    y: float
    z: float
    neighbours: List[int]
    index: int
