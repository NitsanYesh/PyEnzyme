from strenum import StrEnum


class SBOTerm(StrEnum):
    """String enumeration used to assign ontologies derived from SBOTerms."""

    # Chemical reactions
    BIOCHEMICAL_REACTION = "SBO:0000176"
    ACID_BASE_REACTION = "SBO:0000208"
    CONFORMATIONAL_TRANSITION = "SBO:0000181"
    CONVERSION = "SBO:0000182"
    DEGRADATION = "SBO:0000179"
    DISSOCIATION = "SBO:0000180"
    IONISATION = "SBO:0000209"
    ISOMERISATION = "SBO:0000377"
    NON_COVALENT_BINDING = "SBO:0000177"
    REDOX_REACTION = "SBO:0000200"
    SPONTANEOUS_REACTION = "SBO:0000672"

    # Chemical entities
    PROTEIN = "SBO:0000252"
    GENE = "SBO:0000251"
    SMALL_MOLECULE = "SBO:0000247"
    ION = "SBO:0000327"
    RADICAL = "SBO:0000328"

    # Chemical relations
    INTERACTOR = "SBO:0000336"
    SUBSTRATE = "SBO:0000015"
    PRODUCT = "SBO:0000011"
    CATALYST = "SBO:0000013"
    INHIBITOR = "SBO:0000020"
    ESSENTIAL_ACTIVATOR = "SBO:0000461"
    NON_ESSENTIAL_ACTIVATOR = "SBO:0000462"
    POTENTIATOR = "SBO:0000021"

    # Kinetic models
    MICHAELIS_MENTEN = "SBO:0000028"

    # Kinetic parameters
    K_CAT = "SBO:0000025"
    K_M = "SBO:0000027"


class DataTypes(StrEnum):
    """String enumeration used to assign replicate type ontologies"""

    CONCENTRATION = "conc"
    ABSORPTION = "abs"
    FEED = "feed"
    BIOMASS = "biomass"


class EnzymeMLPart(StrEnum):
    """Mapping to identify where entities are stored in the EnzymeML model."""

    # Chemical entities
    PROTEIN = "protein_dict"
    SMALL_MOLECULE = "reactant_dict"
    ION = "reactant_dict"
    RADICAL = "reactant_dict"

    BLANK = "blank"

    @classmethod
    def partFromSBOTerm(cls, sbo_term: str) -> str:
        sbo_term = SBOTerm(sbo_term).name
        return getattr(cls, sbo_term).value

    @classmethod
    def entityFromSBOTerm(cls, sbo_term: str) -> str:
        sbo_term = SBOTerm(sbo_term).name
        return getattr(cls, sbo_term).name
