from ..typing import NamedTuple, List, Optional


class Sequence(NamedTuple):
    """Sequence which contains the index of the monomer in the monomers list, the information if it is inverted and the number of repetitions."""

    index: List[int]
    inverted: List[bool]
    reps: List[int]


class Atom:
    """A class which contains the information about the atom."""

    ff_identifier: str
    """The force field identifier of the atom."""
    element: str
    """The element of the atom."""
    x: float
    """The x coordinate of the atom."""
    y: float
    """The y coordinate of the atom."""
    z: float
    """The z coordinate of the atom."""
    neighbours: List[int]
    """The indices of the neighbouring atoms."""
    index: int
    """The index of the atom in the monomer."""

    def __init__(
        self,
        ff_identifier: Optional[str],
        element: str,
        x: float,
        y: float,
        z: float,
        neighbours: List[int],
        index: int,
    ):
        self.ff_identifier = ff_identifier
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.neighbours = neighbours
        self.index = index
