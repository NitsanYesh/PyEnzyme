'''
File: enzymemldocument.py
Project: core
Author: Jan Range
License: BSD-2 clause
-----
Last Modified: Thursday July 15th 2021 1:00:05 am
Modified By: Jan Range (<jan.range@simtech.uni-stuttgart.de>)
-----
Copyright (c) 2021 Institute of Biochemistry and Technical Biochemistry Stuttgart
'''

import os
import re
import io
import json

from pydantic import Field, validator, PositiveInt, validate_arguments
from typing import TYPE_CHECKING, Optional, Union
from dataclasses import dataclass
from pathlib import Path

from pyenzyme.enzymeml.core.enzymemlbase import EnzymeMLBase
from pyenzyme.enzymeml.core.abstract_classes import AbstractSpecies

from pyenzyme.enzymeml.core.reactant import Reactant
from pyenzyme.enzymeml.core.creator import Creator
from pyenzyme.enzymeml.core.protein import Protein
from pyenzyme.enzymeml.core.vessel import Vessel
from pyenzyme.enzymeml.core.unitdef import UnitDef
from pyenzyme.enzymeml.core.measurement import Measurement
from pyenzyme.enzymeml.core.measurementData import MeasurementData
from pyenzyme.enzymeml.core.enzymereaction import EnzymeReaction
from pyenzyme.enzymeml.tools.unitcreator import UnitCreator
from pyenzyme.enzymeml.tools.enzymemlwriter import EnzymeMLWriter
from pyenzyme.enzymeml.databases.dataverse import toDataverseJSON
from pyenzyme.enzymeml.databases.dataverse import uploadToDataverse

import pyenzyme.enzymeml.tools.enzymemlreader as reader

from pyenzyme.enzymeml.core.ontology import EnzymeMLPart, SBOTerm
from pyenzyme.enzymeml.core.exceptions import SpeciesNotFoundError, IdentifierNameError
from pyenzyme.enzymeml.core.utils import (
    type_checking,
    deprecated_getter
)

from texttable import Texttable


if TYPE_CHECKING:  # pragma: no cover
    static_check_init_args = dataclass
else:
    static_check_init_args = type_checking


@static_check_init_args
class EnzymeMLDocument(EnzymeMLBase):

    name: str = Field(
        description="Title of the EnzymeML Document.",
        required=True
    )

    level: int = Field(
        default=3,
        description="SBML evel of the EnzymeML XML.",
        required=True,
        inclusiveMinimum=1,
        exclusiveMaximum=3
    )

    version: str = Field(
        default=2,
        description="SBML version of the EnzymeML XML.",
        required=True
    )

    pubmed_id: Optional[str] = Field(
        default=None,
        description="Pubmed ID reference.",
        required=False
    )

    url: Optional[str] = Field(
        default=None,
        description="Arbitrary type of URL that is related to the EnzymeML document.",
        required=False,
    )

    doi: Optional[str] = Field(
        default=None,
        description="Digital Object Identifier of the referenced publication or the EnzymeML document.",
        regex=r"/^10.\d{4,9}/[-._;()/:A-Z0-9]+$/i"
    )

    created: Optional[str] = Field(
        default=None,
        description="Date the EnzymeML document was created.",
        required=False
    )

    modified: Optional[str] = Field(
        default=None,
        description="Date the EnzymeML document was modified.",
        required=False
    )

    vessel: Vessel = Field(
        description="The vessel in which the experiment was performed.",
        required=True
    )

    creator_dict: dict[str, Creator] = Field(
        alias="creators",
        default_factory=dict,
        description="Dictionary of the authors of the EnzymeML document.",
        required=True
    )

    protein_dict: dict[str, Protein] = Field(
        alias="proteins",
        default_factory=dict,
        description="Dictionary mapping from protein IDs to protein describing objects.",
        required=False
    )

    reactant_dict: dict[str, Reactant] = Field(
        alias="reactants",
        default_factory=dict,
        description="Dictionary mapping from reactant IDs to reactant describing objects.",
        required=False
    )

    reaction_dict: dict[str, EnzymeReaction] = Field(
        alias="reactions",
        default_factory=dict,
        description="Dictionary mapping from reaction IDs to reaction describing objects.",
        required=False
    )

    unit_dict: dict[str, UnitDef] = Field(
        alias="units",
        default_factory=dict,
        description="Dictionary mapping from unit IDs to unit describing objects.",
        required=False
    )

    measurement_dict: dict[str, Measurement] = Field(
        alias="measurements",
        default_factory=dict,
        description="Dictionary mapping from measurement IDs to measurement describing objects.",
        required=False
    )

    file_dict: dict[str, dict] = Field(
        files="files",
        default_factory=dict,
        description="Dictionary mapping from protein IDs to protein describing objects.",
        required=False
    )

    # ! Validators
    @validator("pubmed_id")
    def add_identifier(cls, pubmed_id: str):
        """Adds an identifiers.org link in front of the pubmed ID if not given"""

        if pubmed_id.startswith("https://identifiers.org/pubmed:"):
            return pubmed_id
        else:
            return "https://identifiers.org/pubmed:" + pubmed_id

    # ! Imports and exports
    @staticmethod
    def fromFile(path: Path):
        """Initializes an EnzymeMLDocument from an OMEX container."

        Args:
            path (Path): Path to the OMEX container.

        Returns:
            EnzymeMLDocument: The intialized EnzymeML document.
        """

        return reader.EnzymeMLReader().readFromFile(path)

    def toFile(self, path: Path, verbose: PositiveInt = 1):
        """Saves an EnzymeML document to an OMEX container at the specified path

        Args:
            path (Path): Path where the document should be saved.
            verbose (PositiveInt, optional): Level of verbosity, in order to print a message and the resulting path. Defaults to 1.
        """

        EnzymeMLWriter().toFile(self, path, verbose=verbose)

    def toXMLString(self):
        """Generates an EnzymeML XML string"""

        return EnzymeMLWriter().toXMLString(self)

    @validate_arguments
    def uploadToDataverse(
        self,
        base_url: str,
        API_Token: str,
        dataverse_name: str
    ):
        """Uploads an EnzymeML document to a Dataverse installation of choice.

        Args:
            base_url (str): URL to a Dataverse installation
            API_Token (str): API Token given from your Dataverse installation for authentication.
            dataverse_name (str): Name of the dataverse to upload the EnzymeML document. You can find the name in the link of your dataverse (e.g. https://dataverse.installation/dataverse/{dataverseName})

        Raises:
            AttributeError: Raised when neither a filename nor an EnzymeMLDocument object was provided.
            ValidationError: Raised when the validation fails.
        """
        uploadToDataverse(
            base_url=base_url,
            API_Token=API_Token,
            dataverse_name=dataverse_name,
            enzmldoc=self
        )

    def toDataverseJSON(self) -> str:
        """Generates a Dataverse compatible JSON representation of this EnzymeML document.

        Returns:
            String: JSON string representation of this EnzymeML document.
        """

        return json.dumps(toDataverseJSON(self), indent=4)

    # ! Utility methods
    @ staticmethod
    def _generateID(prefix: str, dictionary: dict) -> str:
        """Generates IDs complying to the [s|p|r|m|u|c]?[digit]+ schema.

        Args:
            prefix (str): Character denoting the type of species (p: Protein, s: Reactant, u: UnitDef, r: EnzymeReaction, m: Measurement, c: concentration).
            dictionary (dict): The dictionary from which the ID is generated and used to determine the number.

        Returns:
            str: Unique internal identifier.
        """

        if dictionary.keys():
            # fetch all keys and sort them
            number = max(list(dictionary.keys()), key=lambda id: int(id[1::]))
            return prefix + str(number + 1)

        return prefix + str(0)

    def validate(self) -> None:
        # TODO rework validation
        raise NotImplementedError(
            "Function not implemented yet."
        )

    def __str__(self) -> str:
        """
        Magic function return pretty string describing the object.

        Returns:
            string: Beautified summarization of object
        """

        fin_string: list[str]

        def generate_lines(dictionary: dict) -> None:
            """Breaks up a dictionary and generates a human readible line."""
            for element_id, element in dictionary.items():
                fin_string.append(
                    f"\tID: {element_id} \t Name: {element.name}")

        fin_string = ['>>> Units']
        generate_lines(self.unit_dict)

        fin_string.append('>>> Reactants')
        generate_lines(self.reactant_dict)

        fin_string.append('>>> Proteins')
        generate_lines(self.protein_dict)

        fin_string.append('>>> Reactions')
        generate_lines(self.reaction_dict)

        fin_string.append('>>> Measurements')
        fin_string.append(self.printMeasurements())

        return "\n".join(fin_string)

    def printMeasurements(self) -> str:
        """Prints all measurements as a human readable table"""

        table = Texttable()
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["l", "l", "l", "l"])

        # Initialize rows
        rows = [["ID", "Species", "Conc", "Unit"]]

        # Generate and append rows
        for measurement_id, measurement in self.measurement_dict.items():

            speciesDict = measurement.species_dict
            proteins = speciesDict['proteins']
            reactants = speciesDict['reactants']

            # succesively add rows with schema
            # [ measID, speciesID, initConc, unit ]

            for species_id, species in {**proteins, **reactants}.items():
                rows.append(
                    [
                        measurement_id,
                        species_id,
                        str(species.init_conc),
                        species.unit
                    ]
                )

        # Add empty row for better readablity
        rows.append([" "] * 4)
        table.add_rows(rows)

        return f"\n{table.draw()}\n"

    # ! Add methods
    @validate_arguments
    def addReactant(self, reactant: Reactant, use_parser: bool = True) -> str:
        """Adds a Reactant object to the EnzymeML document.

        Args:
            reactant (str): Unique internal identifier of the reactant.
            use_parser (bool, optional): Whether to user the unit parser or not. Defaults to True.

        Returns:
            str: Unique internal identifier of the reactant.
        """

        return self._addSpecies(
            species=reactant,
            prefix="s",
            dictionary=self.reactant_dict,
            use_parser=use_parser,
        )

    @validate_arguments
    def addProtein(self, protein: Protein, use_parser: bool = True) -> str:
        """Adds a Protein object to the EnzymeML document.

        Args:
            protein (str): Unique internal identifier of the reactant.
            use_parser (bool, optional): Whether to user the unit parser or not. Defaults to True.

        Returns:
            str: Unique internal identifier of the reactant.
        """

        return self._addSpecies(
            species=protein,
            prefix="p",
            dictionary=self.protein_dict,
            use_parser=use_parser
        )

    def _addSpecies(
        self,
        species: AbstractSpecies,
        prefix: str,
        dictionary: dict,
        use_parser: bool = True
    ) -> str:
        """Helper function to add any specific species to the EnzymeML document.

        Args:
            species (AbstractSpecies): Species that is about to be added to the EnzymeML document.
            prefix (str): Character that is used to generate a unique internal identifier.
            dictionary (dict): The dictionary where the species will be added to.
            use_parser (bool, optional): Whether to user the unit parser or not. Defaults to True.

        Returns:
            str: The internal identifier of the species.
        """

        # Generate ID
        species.id = self._generateID(
            prefix=prefix, dictionary=dictionary
        )

        # Update unit to UnitDefID
        if species.unit and use_parser:
            species._unit_id = self._convertToUnitDef(species.unit)
        elif species.unit and use_parser is False:
            species._unit_id = species.unit
            species.unit = self.getUnitString(species._unit_id)

        # Add species to dictionary
        dictionary[species.id] = species

        return species.id

    @validate_arguments
    def addReaction(self, reaction: EnzymeReaction, use_parser=True) -> str:
        """
        Adds EnzymeReaction object to EnzymeMLDocument object.
        Automatically assigns ID and converts units.

        Args:
            reaction (EnzymeReaction): Object describing reaction
            use_parser (bool, optional): If set True, will use
                                         internal unit parser.
                                         Defaults to True.

        Returns:
            string: Internal identifier for the reaction.
            Use it for other objects!
        """

        # Generate ID
        reaction.id = self._generateID("r", self.reaction_dict)

        if use_parser:
            # Reset temperature for SBML compliance to Kelvin
            reaction.temperature = (
                reaction.temperature + 273.15
                if re.match(r"^K|kelvin", reaction.temperature_unit)
                else reaction.temperature
            )

            # Generate internal ID for the unit
            reaction._temperature_unit_id = self._convertToUnitDef(
                reaction.temperature_unit
            )
        else:
            # Set the temperature unit to the actual string
            reaction._temperature_unit_id = reaction.temperature
            reaction.temperature_unit = self.getUnitString(
                reaction.temperature_unit
            )

        # Set model units
        if hasattr(reaction, '_EnzymeReaction_model'):
            # TODO implement model unit conversion
            model = reaction.getModel()

        # Finally add the reaction to the document
        self.reaction_dict[reaction.id] = reaction

        return reaction.id

    def addFile(
        self,
        filepath=None,
        fileHandle=None,
        description="Undefined"
    ) -> str:
        """Adds any arbitrary file to the document. Please note, that if a filepath is given, any fileHandle will be ignored.

        Args:
            filepath (str, optional): Path to the file that is added to the document. Defaults to None.
            fileHandle (io.BufferedReader, optional): File handle that will be read to a bytes string. Defaults to None.

        Returns:
            str: Internal identifier for the file.
        """

        # Generate a unique identifier for the file
        file_id = self._generateID("f", self.file_dict)

        if filepath:
            # Open file handle
            fileHandle = open(filepath, "rb")
        elif filepath is None and fileHandle is None:
            raise ValueError(
                "Please specify either a file path or a file handle"
            )

        # Finally, add the file and close the handler
        self.file_dict[file_id] = {
            "name": os.path.basename(fileHandle.name),
            "content": fileHandle.read(),
            "description": description
        }

        fileHandle.close()

        return file_id

    @validate_arguments
    def addMeasurement(self, measurement: Measurement) -> str:
        """Adds a measurement to an EnzymeMLDocument and validates consistency with already defined elements of the documentself.

        Args:
            measurement (Measurement): Collection of data and initial concentrations per reaction

        Returns:
            measurement_id (String): Assigned measurement identifier.
        """

        # Check consistency
        self._checkMeasurementConsistency(measurement)

        # Convert all measurement units to UnitDefs
        self._convertMeasurementUnits(measurement)

        # Finally generate the ID and add it to the dictionary
        measurement_id = self._generateID(
            prefix="m", dictionary=self.measurement_dict
        )
        measurement.id = measurement_id

        self.measurement_dict[measurement_id] = measurement

        return measurement_id

    def _convertMeasurementUnits(self, measurement: Measurement) -> None:
        """Converts string SI units to UnitDef objects and IDs

        Args:
            measurement (Measurement): Object defining a measurement
        """
        species_dict = measurement.species_dict
        measurement.global_time_unit = self._convertToUnitDef(
            measurement.global_time_unit
        )

        def update_unit(measurement_data: MeasurementData) -> None:
            """Helper function to update units"""
            unit = measurement_data.unit
            measurement_data.unit = self._convertToUnitDef(unit)
            self._convertReplicateUnits(measurement_data)

        # Perform update
        map(update_unit, species_dict["proteins"].values())
        map(update_unit, species_dict["reactants"].values())

    def _convertReplicateUnits(self, measurement_data: MeasurementData) -> None:
        """Converts replicate unit strings to unit definitions.

        Args:
            measurement_data (MeasurementData): Object holding measurement data for a species
        """
        for replicate in measurement_data.replicates:

            # Convert unit
            time_unit = self._convertToUnitDef(replicate.time_unit)
            data_unit = self._convertToUnitDef(replicate.data_unit)

            # Assign unit IDs
            replicate.data_unit = data_unit
            replicate.time_unit = time_unit

    def _checkMeasurementConsistency(self, measurement: Measurement) -> None:
        """Checks if the used species in the measurement are consistent with the EnzymeML document.

        Args:
            measurement (MeasurementData): Objech holding measurement data for a species.
        """

        map(self._checkSpecies, measurement.species_dict["reactants"])
        map(self._checkSpecies, measurement.species_dict["proteins"])

    def _checkSpecies(self, species_id):
        """Checks if a species is defined in the EnzymeML document.

        Args:
            species_id (str): Unique identifier of the species.

        Raises:
            SpeciesNotFoundError: Raised when a species is not defined in the EnzymeML document.
        """

        all_species = {**self.reactant_dict, **self.protein_dict}

        if species_id not in all_species.keys():

            # Retrieve species for ontology
            species = self._getSpecies(
                id=species_id,
                dictionary=all_species,
                element_type="Proteins/Reactants"
            )

            # Use the EnzymeMLPart Enum to derive the correct place
            sbo_term = SBOTerm(species.__dict__["ontology"]).name
            enzymeml_part = EnzymeMLPart.fromSBOTerm(sbo_term)

            # Raise an error if the species is nowhere present
            raise SpeciesNotFoundError(
                species_id=species_id,
                enzymeml_part=enzymeml_part
            )

    def _convertToUnitDef(self, unit: Optional[str]) -> str:
        """Reads an SI unit string and converts it into a EnzymeML compatible UnitDef

        Args:
            unit (str): String representing the SI unit.

        Returns:
            str: Unique identifier of the UnitDef.
        """
        if unit is None:
            raise TypeError("No unit given.")
        elif unit in self.unit_dict.keys():
            return unit

        return UnitCreator().getUnit(unit, self)

    # ! Getter methods

    def getUnitString(self, unit_id: str) -> str:
        """Return the unit name corresponding to the given unit ID.

        Args:
            unit_id (str): Unique internal ID of the unit.

        Raises:
            SpeciesNotFoundError: Raised when the requested unit is not found.

        Returns:
            str: String representation of the unit.
        """

        try:
            return self.unit_dict[unit_id].name
        except KeyError:
            raise SpeciesNotFoundError(
                species_id=unit_id, enzymeml_part="Units"
            )

    def getUnitDef(self, id: str, by_id: bool = True):
        """Returns the unit associated with the given ID.

        Args:
            id (str): Unique internal ID of the unit.
            by_id (bool, optional): Whether the unit is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested unit is not found.

        Returns:
            UnitDef: The corresponding unit object.
        """

        self._getSpecies(
            id=id,
            dictionary=self.unit_dict,
            element_type="Units",
            by_id=by_id
        )

    def getReaction(self, id: str, by_id: bool = True):
        """Returns the reaction associated with the given ID.

        Args:
            id (str): Unique internal ID of the reaction.
            by_id (bool, optional): Whether the reaction is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested reaction is not found.

        Returns:
            EnzymeReaction: The corresponding reaction object.
        """

        return self._getSpecies(
            id=id,
            dictionary=self.reaction_dict,
            element_type="EnzymeReaction",
            by_id=by_id
        )

    def getMeasurement(self, id: str, by_id: bool = True):
        """Returns the measurement associated with the given ID.

        Args:
            id (str): Unique internal ID of the measurement.
            by_id (bool, optional): Whether the measurement is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested measurement is not found.

        Returns:
            Measurement: The corresponding measurement object.
        """

        return self._getSpecies(
            id=id,
            dictionary=self.measurement_dict,
            element_type="Measurement",
            by_id=by_id
        )

    def getReactant(self, id: str, by_id=True):
        """Returns the reactant associated with the given ID.

        Args:
            id (str): Unique internal ID of the reactant.
            by_id (bool, optional): Whether the reactant is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested reactant is not found.

        Returns:
            Reactant: The corresponding reactant object.
        """

        return self._getSpecies(
            id=id,
            dictionary=self.reactant_dict,
            element_type="Reactant",
            by_id=by_id
        )

    def getProtein(self, id: str, by_id: bool = True):
        """Returns the protein associated with the given ID.

        Args:
            id (str): Unique internal ID of the protein.
            by_id (bool, optional): Whether the protein is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested protein is not found.

        Returns:
            Protein: The corresponding protein object.
        """

        return self._getSpecies(
            id=id,
            dictionary=self.protein_dict,
            element_type="Protein",
            by_id=by_id
        )

    def getFile(self, id: str, by_id: bool = True):
        """Returns the file associated with the given ID.

        Args:
            id (str): Unique internal ID of the file.
            by_id (bool, optional): Whether the file is retrieved via ID or name. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested file is not found.

        Returns:
            dict[str, dict]: The corresponding file object.
        """

        return self._getSpecies(
            id=id,
            dictionary=self.file_dict,
            element_type="File",
            by_id=by_id
        )

    def _getSpecies(
        self,
        id: str,
        dictionary: dict,
        element_type: str,
        by_id: bool = True
    ) -> Union[
            AbstractSpecies,
            EnzymeReaction,
            Measurement,
            UnitDef,
            dict
    ]:
        """Helper function to retrieve any kind of species from the EnzymeML document.

        Args:
            id (str): Unique internal ID.
            dictionary (dict): Dictionary that stores all objects.
            element_type (str): Type of object that is in the dictionary.
            by_id (bool, optional): [description]. Defaults to True.

        Raises:
            SpeciesNotFoundError: Raised when the requested species is not found.

        Returns:
            Union[ AbstractSpecies, EnzymeReaction, Measurement ]: The requested object
        """

        # Fix the searched attribute
        searched_attribute = "id" if by_id else "name"

        try:
            # Filter the dict for the desired species
            return next(filter(
                lambda obj: obj.__dict__[searched_attribute] == id,
                dictionary.values()
            ))
        except StopIteration:
            # When the generator is empty, raise error
            raise SpeciesNotFoundError(
                species_id=id, enzymeml_part=element_type
            )

    def getReactantList(self) -> list[Reactant]:
        """Returns a list of all reactants in the EnzymeML document."

        Returns:
            list[Reactant]: List of all reactants in the EnzymeML document.
        """
        return self._getSpeciesList(self.reactant_dict)

    def getProteinList(self) -> list[Protein]:
        """Returns a list of all proteins in the EnzymeML document."

        Returns:
            list[Protein]: List of all proteins in the EnzymeML document.
        """
        return self._getSpeciesList(self.protein_dict)

    def getReactionList(self) -> list[EnzymeReaction]:
        """Returns a list of all reactions in the EnzymeML document."

        Returns:
            list[EnzymeReaction]: List of all reactions in the EnzymeML document.
        """
        return self._getSpeciesList(self.reaction_dict)

    def getFilesList(self):
        """Returns a list of all files in the EnzymeML document."

        Returns:
            list[dict]: List of all files in the EnzymeML document.
        """
        return self._getSpeciesList(self.file_dict)

    @ staticmethod
    def _getSpeciesList(dictionary: dict) -> list:
        """Helper function to retrieve lists of dicitonary objects

        Args:
            dictionary (dict): Dictionary of corresponding elements

        Returns:
            list: Returns all values in the dictionary
        """
        return list(dictionary.values())

    @deprecated_getter("doi")
    def getDoi(self) -> Optional[str]:
        return self.doi

    @deprecated_getter("pubmed_id")
    def getPubmedID(self) -> Optional[str]:
        return self.pubmed_id

    @deprecated_getter("url")
    def getUrl(self) -> Optional[str]:
        return self.url

    @deprecated_getter("created")
    def get_created(self):
        return self.created

    @deprecated_getter("modified")
    def getModified(self):
        return self.modified

    @deprecated_getter("creators")
    def getCreator(self):
        return self.creator_dict

    @deprecated_getter("vessel")
    def getVessel(self):
        return self.vessel

    @deprecated_getter("name")
    def getName(self):
        return self.name

    @deprecated_getter("level")
    def getLevel(self):
        return self.level

    @deprecated_getter("version")
    def getVersion(self):
        return self.version

    @deprecated_getter("protein_dict")
    def getProteinDict(self):
        return self.protein_dict

    @deprecated_getter("reactant_dict")
    def getReactantDict(self):
        return self.reactant_dict

    @deprecated_getter("reaction_dict")
    def getReactionDict(self):
        return self.reaction_dict

    @deprecated_getter("measurement_dict")
    def getMeasurementDict(self):
        return self.measurement_dict

    @deprecated_getter("unit_dict")
    def getUnitDict(self):
        return self.unit_dict

    @deprecated_getter("file_dict")
    def getFileDict(self):
        return self.file_dict
