import numpy as np
import pandas as pd
import re

from typing import Any, List, Dict, Tuple
from math import isnan
from collections import OrderedDict

from pyenzyme.enzymeml.core.ontology import DataTypes
from pyenzyme.enzymeml.core.creator import Creator
from pyenzyme.enzymeml.core.vessel import Vessel
from pyenzyme.enzymeml.core.protein import Protein
from pyenzyme.enzymeml.core.reactant import Reactant
from pyenzyme.enzymeml.core.enzymereaction import EnzymeReaction
from pyenzyme.enzymeml.core.measurement import Measurement
from pyenzyme.enzymeml.core.replicate import Replicate
from pyenzyme.enzymeml.tools.templatereader import get_instances, add_instances, merge_protein_modifier, parse_reaction_element


def read_96well_template(path: str, enzmldoc):

    general_info = pd.read_excel(
        path, sheet_name="General Information", skiprows=1)

    params = dict(
        name=general_info.iloc[0, 1],
        created=str(general_info.iloc[1, 1]),
        doi=None,
        pubmedid=general_info.iloc[3, 1],
        url=general_info.iloc[4, 1],
    )

    enzmldoc = enzmldoc(**params)

    # User information
    user_infos = pd.read_excel(
        path, sheet_name="General Information", skiprows=9, nrows=8).dropna()

    for record in user_infos.to_dict(orient="records"):
        enzmldoc.addCreator(
            Creator(
                family_name=record["Family Name"],
                given_name=record["Given Name"],
                mail=record["Mail"],
            )
        )

    # Vessel (96 well plate)
    vessel_info = pd.read_excel(
        path, sheet_name="General Information", skiprows=18, nrows=2).dropna()

    vessel_info = vessel_info.to_dict(orient="records")[0]

    vessel_id = enzmldoc.addVessel(
        Vessel(
            id=vessel_info["ID"],
            name=vessel_info["Name"],
            volume=vessel_info["Volume value"],
            unit=vessel_info["Volume unit"]
        )
    )
    vessel = {"vessel_id": vessel_id}

    # Proteins
    proteins = pd.read_excel(path, sheet_name="Proteins", skiprows=2)
    instances = get_instances(proteins, Protein, enzmldoc)

    for instance in instances:
        enzmldoc.addProtein(Protein(**instance | vessel))

    # Reactants
    reactants = pd.read_excel(
        path, sheet_name="Chemicals", skiprows=2, usecols="A:E")
    instances = get_instances(reactants, Reactant, enzmldoc)

    for instance in instances:
        enzmldoc.addReactant(Reactant(**instance | vessel))

    # Reactions
    reactions = pd.read_excel(
        path, sheet_name="Reactions", skiprows=2, usecols="A:J")

    # Merge proteins and modifiers
    nu_mods = [
        merge_protein_modifier(protein, modifier)
        for modifier, protein in zip(
            reactions.Modifiers.values.tolist(), reactions.Proteins.values.tolist()
        )
    ]

    # Replace merged modifiers with modifier tag
    reactions.Modifiers = nu_mods

    instances = get_instances(reactions, EnzymeReaction, enzmldoc)

    for instance in instances:

        # Get Educts, Products and Modifiers to add to the reaction
        educts = parse_reaction_element(instance.get("educts"))
        products = parse_reaction_element(instance.get("products"))
        modifiers = parse_reaction_element(instance.get("modifiers"))

        instance.pop("educts")
        instance.pop("products")
        instance.pop("modifiers")

        # Instantiate Reaction
        reaction = EnzymeReaction(**instance)

        add_instances(reaction.addModifier, modifiers, enzmldoc)
        add_instances(reaction.addEduct, educts, enzmldoc)
        add_instances(reaction.addProduct, products, enzmldoc)

        enzmldoc.addReaction(reaction)

    # Set initial conditions of measurements
    enzmldoc = generate_measurements(path, enzmldoc)

    for measurement in enzmldoc.measurement_dict.values():

        # get initial concentrations of reactants
        for reactant_id in enzmldoc.reactant_dict:
            reactant_name = enzmldoc.getReactant(reactant_id).name
            unit = get_species_unit(path, reactant_name)
            init_concs = extract_initial_conditions(
                path, reactant_name)

            measurement.addData(
                init_conc=init_concs[measurement.name],
                unit=unit,
                reactant_id=reactant_id
            )

        # get initial concentrations of proteins
        for protein_id in enzmldoc.protein_dict:
            protein_name = enzmldoc.getProtein(protein_id).name
            unit = get_species_unit(path, protein_name)
            init_concs = extract_initial_conditions(
                path, protein_name)

            measurement.addData(
                init_conc=init_concs[measurement.name],
                unit=unit,
                protein_id=protein_id
            )

    # Add measurement data
    type_mapping = {
        "Concentration": DataTypes.CONCENTRATION,
        "Absorption": DataTypes.ABSORPTION,
        "Conversion [%]": DataTypes.CONVERSION,
        "Peak Area": DataTypes.PEAK_AREA,
        "Total concentration after addition": DataTypes.CONCENTRATION,
    }

    data_info = pd.read_excel(path, skiprows=2, sheet_name="Data", nrows=0)
    measured_reactant = data_info.columns[2]
    data_type = data_info.columns[5]
    time_unit = data_info.columns[8]

    measured_reactant_id = [reactant.id for reactant in enzmldoc.reactant_dict.values(
    ) if reactant.name == measured_reactant][0]

    measured_data = pd.read_excel(path, skiprows=3, sheet_name="Data")
    time = measured_data.pop("Time").values.tolist()
    data_dict = measured_data.to_dict(orient="list")

    for measurement in enzmldoc.measurement_dict.values():

        measurement.addReplicates(
            Replicate(
                id=measurement.name,
                species_id=measured_reactant_id,
                time_unit=time_unit,
                time=time,
                data=data_dict[measurement.name],
                is_calculated=False,
                data_type=type_mapping[data_type],
                data_unit="dimensionless"
            ),
            enzmldoc
        )
    return enzmldoc


def validate_plate_layout_homogeneity(path: str, enzmldoc):
    """Validates that for all species the initial concentration is set
    for all well positions."""

    proteins = [protein.name for protein in enzmldoc.protein_dict.values()]
    reactants = [reactant.name for reactant in enzmldoc.reactant_dict.values()]

    species = proteins + reactants

    well_ids = []
    species_well_ids = []
    for spec in species:
        spec_well_ids = extract_initial_conditions(path, spec).keys()
        species_well_ids.append(spec_well_ids)
        well_ids = well_ids + list(spec_well_ids)

    unique_well_ids = set(well_ids)

    for spec, spec_well_ids in zip(species, species_well_ids):
        if set(spec_well_ids) != unique_well_ids:
            raise ValueError(
                f"Initial concentration for plate postion(s)"
                f"{list(set(spec_well_ids) ^ set(unique_well_ids))} of {spec} not specified."
            )


def extract_initial_conditions(path: str, sheet_name: str) -> OrderedDict[str, float]:
    """Extracts initial concentration values for each cell, representing a 
    position on a 96 well plate."""

    plate = pd.read_excel(path, sheet_name=sheet_name, skiprows=3,
                          nrows=8, usecols='A:M', header=0, index_col=0)

    plate_dict = plate.T.to_dict(orient='dict')

    # flatten nested dict
    flat_plate_dict = flatten_dict(plate_dict)

    # remove nan values
    clean_dict = OrderedDict({key: value for key, value in flat_plate_dict.items(
    ) if not isnan(flat_plate_dict[key])})

    return clean_dict


def flatten_dict(nested_dict: Dict[str, Dict]) -> dict:
    """Flattens nested dict. Concatenates inner and outer keys."""

    flat_dict = OrderedDict()
    for outer_k, outer_v in nested_dict.items():
        for inner_k, inner_v in outer_v.items():
            flat_dict.update({outer_k.strip()+str(inner_k): inner_v})

    return flat_dict


def get_timecourse_data(path: str) -> Tuple[List[float], Dict[str, List[float]]]:

    df = pd.read_excel(path, skiprows=3, sheet_name="Data")
    time = df.pop("Time").values.tolist()
    data_dict = df.to_dict(orient="list")

    return time, data_dict


def get_species_unit(path: str, sheet_name: str) -> str:
    """Extracts the unit of the reactant / protein from the template."""

    unit = pd.read_excel(path, sheet_name=sheet_name,
                         skiprows=2, nrows=0, usecols="D")
    return unit.columns[0]


def generate_measurements(path: str, enzmldoc):
    """Generates a measurement for each well position in the template."""

    validate_plate_layout_homogeneity(path, enzmldoc)

    reactant_name = next(iter(enzmldoc.reactant_dict.values())).name

    well_positions = extract_initial_conditions(path, reactant_name).keys()
    for well in well_positions:
        enzmldoc.addMeasurement(Measurement(name=well))

    return enzmldoc
