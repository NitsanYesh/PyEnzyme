from pydantic import Field, BaseModel
from typing import Optional, List

from pyenzyme.enzymeml.core.enzymemlbase import EnzymeMLBase

# We can design this in different ways:
# Eighter calculated data is provided in the form of concentrations 
# (e.g. the stock solution consisted of 0.1 mol/l s0 and 50 mmol/l s1),
#
# or we document "labratory raw-data". Then pyenzyme calculates the concentration based on mass, molar_mass (, and density).
# (e.g. 30 mg of s1 were solved in 5 mL of s0) is provided in the form



class StockSolutionElement(BaseModel):
    """Describes an element of an reaction mixture."""

    species_id: str = Field( # Can be reactant, or protein.
        ...,
        description="Internal identifier to either a protein or reactant defined in the EnzymeMLDocument.",
    )

    # Depending on the state of matter, eighter mass or volume needs to be provided.
    mass: Optional[float] = Field(
        None,
        description="Mass of the reactant in the stock solution."
    )

    mass_unit: Optional[str] = Field(
        None,
        description="Unit of the reactants mass."
    )

    volume: Optional[float] = Field(
        None,
        description="Volume of the reactant in the stock solution."
    )

    volume_unit: Optional[str] = Field(
        None,
        description="Unit of the reactants volume in the stock solution."
    )


class StockSolution(EnzymeMLBase):

    name: Optional[str] = Field(
        None, description="Name of the stock solution."
    )

    # I guess vessel can be optional, since it only serves as a temporary container.
    vessel_id: Optional[str] = Field(
        description="Identifier of the vessel in which the stock solution was stored.",
        template_alias="Vessel",
    )

    components: List[StockSolutionElement] = Field(
        default_factory=list,
        description="List of reactants containing StockSolutionElement objects."
    )


    def addComponent(
        self,
        species_id: str,
    ) -> None:
        pass


    @classmethod
    def getConcentration(self):
        pass