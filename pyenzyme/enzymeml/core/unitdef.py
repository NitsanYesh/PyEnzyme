'''
File: unitdef.py
Project: core
Author: Jan Range
License: BSD-2 clause
-----
Last Modified: Tuesday June 15th 2021 7:48:57 pm
Modified By: Jan Range (<jan.range@simtech.uni-stuttgart.de>)
-----
Copyright (c) 2021 Institute of Biochemistry and Technical Biochemistry Stuttgart
'''

from pydantic import Field, validator, validate_arguments
from typing import TYPE_CHECKING, Optional
from dataclasses import dataclass

from pyenzyme.enzymeml.core.enzymemlbase import EnzymeMLBase
from pyenzyme.enzymeml.core.utils import (
    type_checking,
    deprecated_getter
)


if TYPE_CHECKING:  # pragma: no cover
    static_check_init_args = dataclass
else:
    static_check_init_args = type_checking


@static_check_init_args
class BaseUnit(EnzymeMLBase):
    """Base unit description including kind, exponent, scale and multiplier"""

    kind: str = Field(
        description="Unit kind used to write SBML.",
        required=True
    )

    exponent: float = Field(
        description="Unit exponent.",
        required=True
    )

    scale: int = Field(
        description="Unit scale.",
        required=True
    )

    multiplier: float = Field(
        description="Unit multiplier.",
        required=True
    )

    # @validator("kind")
    # def check_sbml_unit_enum(cls, kind_int: int):
    #     kind_string: str = libsbml.UnitKind_toString(kind_int)

    #     if "Invalid UnitKind" in kind_string:
    #         raise UnitKindError()


@static_check_init_args
class UnitDef(EnzymeMLBase):

    name: str = Field(
        description="Name of the SI unit.",
        required=True
    )

    id: Optional[str] = Field(
        default=None,
        description="Interal Identifier of the SI unit.",
        required=True
    )

    meta_id: Optional[str] = Field(
        description="Interal meta identifier of the SI unit.",
        required=True
    )

    units: list[BaseUnit] = Field(
        default_factory=list,
        description="List of SI baseunits.",
        required=True
    )

    ontology: Optional[str] = Field(
        default=None,
        description="Ontology of the SI unit.",
        required=False
    )

    # ! Validators
    @validator("id")
    def set_meta_id(cls, id: Optional[str], values: dict):
        """Sets the meta ID when an ID is provided"""

        if id:
            # Set Meta ID with ID
            values["meta_id"] = f"METAID_{id.upper()}"

        return id

    @validator("meta_id")
    def check_meta_id(cls, meta_id: Optional[str], values: dict):
        """Checks if the meta ID provided is following the standard"""

        if values.get("meta_id"):
            # When the ID init already set the meta ID
            return values.get("meta_id")

        return None

    @validate_arguments
    def addBaseUnit(self, kind: str, exponent: float, scale: int, multiplier: float) -> None:
        """Adds a base unit to the units element and sort the units.

        Args:
            kind (str): SBML unit kind string.
            exponent (float): Exponent of the unit.
            scale (float): Scale of the unit.
            multiplier (float): Muliplier of the unit.
        """

        # Create baseunit
        baseunit = BaseUnit(kind=kind, exponent=exponent,
                            scale=scale, multiplier=multiplier)

        # Merge both and sort them via kind
        self.units.append(baseunit)
        self.units = sorted(
            self.units,
            key=lambda unit: unit.kind
        )

    # ! Utilities
    def calculateTransformValue(self, kind: str, scale: int):
        """Calculates the value that is needed to re-scale the given unit to the desired scale.

        Args:
            kind (str): The kind of unit that is used as a reference.
            scale (int): The desired scale.

        Raises:
            ValueError: Raised when the given unit kind is not part of the unitdef.

        Returns:
            float: The value that is needed to re-scale the given unit to the desired scale.
        """

        for base_unit in self.units:
            if base_unit.kind == kind:

                # correction factor used for the case of scale=1
                coorrection_factor = -1 if base_unit.scale == 1 else 0

                return 10 ** (
                    base_unit.exponent * (
                        base_unit.scale - scale + coorrection_factor
                    )
                )

        raise ValueError(
            f"Unit kind of {kind} is not part of the unit definition"
        )

    def _getNewName(self) -> str:
        """Internal function used to derive a units new name. Will be assigned using enzmldoc._convertTounitDef.

        Returns:
            str: The new name of the unit definition.
        """

        # Mapping for abbreviations
        kind_mapping = {
            "mole": "mole",
            "second": "s",
            "liter": "l",
            "litre": "l",
        }

        prefix_mapping = {
            -15: "f",
            -12: "p",
            -9: "n",
            -6: "u",
            -3: "m",
            -2: "c",
            -1: "d",
            1: "",
            3: "k"
        }

        nominator = list(filter(
            lambda base_unit: base_unit.exponent > 0,
            self.units
        ))

        denominator = list(filter(
            lambda base_unit: base_unit.exponent < 0,
            self.units
        ))

        # Create new unit name
        def constructName(base_unit: BaseUnit) -> str:
            return f"{prefix_mapping[base_unit.scale]}{kind_mapping[base_unit.kind]}"

        nominator_string = " ".join([constructName(base_unit)
                                     for base_unit in nominator])
        denominator_string = " ".join([constructName(base_unit)
                                       for base_unit in denominator])

        return " / ".join([nominator_string, denominator_string])

    # ! Getters
    @deprecated_getter("units")
    def getUnits(self):
        return self.units

    @deprecated_getter("name")
    def getName(self):
        return self.name

    @deprecated_getter("id")
    def getId(self):
        return self.id

    @deprecated_getter("meta_id")
    def getMetaid(self):
        return self.meta_id

    @deprecated_getter("ontology")
    def getOntology(self):
        return self.ontology

    def getFootprint(self):
        sorted_units = [base_unit.dict() for base_unit in self.units]
        return sorted(
            sorted_units,
            key=lambda unit: unit["kind"]
        )
