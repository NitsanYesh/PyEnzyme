from pydantic import Field, BaseModel
from typing import Optional, List, Dict

from pyenzyme.enzymeml.core.abstract_classes import AbstractSpecies
from pyenzyme.enzymeml.core.reactant import Reactant
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

    components_dict: Dict[str, StockSolutionElement] = Field(
        default_factory=list,
        description="List of reactants containing StockSolutionElement objects."
    )

    def addComponent(
        self,
        species_id: str,
        mass: Optional[float] = None,
        mass_unit: Optional[str] = None,
        volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
    ) -> None:

        self.components_dict[species_id] = StockSolutionElement(
            species_id=species_id,
            mass=mass,
            mass_unit=mass_unit,
            volume=volume,
            volume_unit=volume_unit
        )



    def getConcentrations(self, reactant_dict: Dict[str, Reactant]) -> Dict[str, float]:

        concentration_in_stock_solution_dict = dict()
        total_volume = self._getTotalVolume()

        for reactant_id, reactant in self.components_dict.items():
            
            # Check if all species have molar mass
            if reactant_dict[reactant_id].molar_mass is not None:
                molar_mass = reactant_dict[reactant_id].molar_mass
            else:
                raise TypeError(f"No molar mass specified for species {reactant_id}.")

            # Calaculate concentration in stock solution of dissolved solids
            if not self.components_dict[reactant_id].mass == None:

                mols = reactant.mass / molar_mass
                concentration = mols / total_volume
                concentration_in_stock_solution_dict[reactant_id] = concentration
            
            # Calculate concentration of fluids
            if not self.components_dict[reactant_id].volume == None:

                # Check if density is provided for all fluids
                if reactant_dict[reactant_id].density is not None:
                    density = reactant_dict[reactant_id].density
                else:
                    raise TypeError(f"No density specified for species {reactant_id}.")

                volumetric_molarity = molar_mass / density # l / mol
                mols = reactant.volume / volumetric_molarity
                concentration = mols / total_volume
                concentration_in_stock_solution_dict[reactant_id] = concentration




        return concentration_in_stock_solution_dict

    def _getTotalVolume(self) -> float:
        """Calculate total volume of stock solution"""

        total_volume = 0.0
        for component in self.components_dict.values():
            if component.volume != None:
                total_volume += component.volume

        return total_volume





if __name__ == "__main__":
    from pyenzyme.enzymeml.core.reactant import Reactant

    reactant_dict = {}

    # Fluid components (water and methanol)
    reactant_dict["s0"] = Reactant.fromChebiID(15377, "v0")
    reactant_dict["s0"].density = 997

    reactant_dict["s1"] = Reactant.fromChebiID(17790, "v0")
    reactant_dict["s1"].density = 792
    
    # Solid components
    reactant_dict["s2"] = Reactant.fromChebiID(27732, "v0")


    stocksolution = StockSolution(name="Substrate stock")
    
    stocksolution.addComponent("s0", volume=0.005, volume_unit="l")
    stocksolution.addComponent("s1", volume=0.0001, volume_unit="l")
    stocksolution.addComponent("s2", mass=0.005, mass_unit="g")

    print(stocksolution.getConcentrations(reactant_dict))

