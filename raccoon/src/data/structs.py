from ..typing import NamedTuple, List


class Sequence(NamedTuple):
    """Sequence which contains the index of the monomer in the monomers list, the information if it is inverted and the number of repetitions."""

    index: List[int]
    inverted: List[bool]
    reps: List[int]


class Atom(NamedTuple):
    """A Namedtuple which contains the information about the atom."""

    name: str
    index: int
    x: float
    y: float
    z: float
    ff_identifier: str
    neighbours: List[int]
