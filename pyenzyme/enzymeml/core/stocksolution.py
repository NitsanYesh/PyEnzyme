from pydantic import Field, BaseModel
from typing import Optional, List, Dict

from pyenzyme.enzymeml.core.abstract_classes import AbstractSpecies
from pyenzyme.enzymeml.core.reactant import Reactant
from pyenzyme.enzymeml.core.enzymemlbase import EnzymeMLBase
from pyenzyme.enzymeml.core.exceptions import SpeciesNotFoundError


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
        None,
        description="Name of the stock solution."
    )

    id: Optional[str] = Field(
        None,
        description="Unique identifier of the stock solution."
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

    concentration_dict: Dict[str, float] = Field(
        None,
        description="Calculated concentrations of each component within the stock solution."
    )

    def addComponent(
        self,
        species_id: str,
        mass: Optional[float] = None,
        mass_unit: Optional[str] = None,
        volume: Optional[float] = None,
        volume_unit: Optional[str] = None,
    ) -> None:

        self.components.append(StockSolutionElement(
            species_id=species_id,
            mass=mass,
            mass_unit=mass_unit,
            volume=volume,
            volume_unit=volume_unit
        ))

    def getConcentrations(self, all_dict: Dict[str, Reactant]) -> Dict[str, float]:
        """ Calculates molar concentration of each component within the stock solution.
        
        Args:
            all_dict (dict): dict, all enzymeML dicts.

        Returns:
             concentration_in_stock_solution_dict (dict): reaction_id: molar concentration pairs.
        """

        # TODO handling of units
        # TODO dict.value should store aoncentration as well as the respective concentration

        concentration_in_stock_solution_dict = dict()
        total_volume = self._getTotalVolume()

        for species in self.components:
            species_id = self._getSpecies(species.species_id).species_id
            
            # Check if all species have molar mass
            if all_dict[species_id].molar_mass is not None:
                molar_mass = all_dict[species_id].molar_mass
            else:
                raise ValueError(f"No molar mass specified for species {species_id}.")

            # Calaculate concentration in stock solution of dissolved solids
            if species.mass is not None:

                mols = species.mass / molar_mass
                concentration = mols / total_volume
                concentration_in_stock_solution_dict[species_id] = concentration
            
            # Calculate concentration of fluids
            if species.volume is not None:

                # Check if density is provided for all fluids
                if all_dict[species_id].density is not None:
                    density = all_dict[species_id].density
                else:
                    raise ValueError(f"No density specified for species {species_id}.")

                volumetric_molarity = molar_mass / density # l / mol
                mols = species.volume / volumetric_molarity
                concentration = mols / total_volume
                concentration_in_stock_solution_dict[species_id] = concentration

        return concentration_in_stock_solution_dict

    def _getTotalVolume(self) -> float:
        """Calculate total volume of stock solution"""

        total_volume = 0
        for component in self.components:
            if component.volume is not None:
                total_volume += component.volume

        return total_volume


    def _getSpecies(self, id: str) -> "StockSolutionElement":
        try:
            return next(filter(lambda element: element.species_id == id, self.components))
        except StopIteration:
            raise SpeciesNotFoundError(species_id=id, enzymeml_part="StockSolution")

    @staticmethod 
    def _calculate_concentration(
        stock_solution: "StockSolution",
        volume: float
        ) -> float:
        pass








if __name__ == "__main__":
    from pyenzyme.enzymeml.core.reactant import Reactant
    from pyenzyme.enzymeml.core.protein import Protein
    from pyenzyme.enzymeml.core.measurement import Measurement


    reactant_dict = {}
    enzyme_dict = {}
    stocksolution_dict = {}

    ###################
    ## Fluid Species ##
    ################### 
      
    # Water
    reactant_dict["s0"] = Reactant.fromChebiID(15377, "v0")
    reactant_dict["s0"].density = 997

    # Ethanol
    reactant_dict["s1"] = Reactant.fromChebiID(16236, "v0")
    reactant_dict["s1"].density = 789
    

    ###################
    ## Solid Species ##
    ###################

    # Coffein
    reactant_dict["s2"] = Reactant.fromChebiID(27732, "v0")

    # Tris
    reactant_dict["s3"] = Reactant.fromChebiID(9754, "v0")

    #############
    ## Protein ##
    #############

    enzyme_dict["p0"] = Protein.fromUniProtID("Q16678", "v0")
    enzyme_dict["p0"].molar_mass = 60860


    ###########################
    ## Define StockSolutions ##
    ###########################  

    # Substrate stock
    stocksolution_dict["sto0"] = StockSolution(name="Disco Mate")
    stocksolution_dict["sto0"].addComponent("s0", volume=0.450, volume_unit="l")
    stocksolution_dict["sto0"].addComponent("s1", volume=0.05, volume_unit="l")
    stocksolution_dict["sto0"].addComponent("s2", mass=0.02, mass_unit="g")

    # Tris-Buffer stock
    stocksolution_dict["sto1"] = StockSolution(name="Hepes Buffer")
    stocksolution_dict["sto1"].addComponent("s0", volume=1, volume_unit="l")
    stocksolution_dict["sto1"].addComponent("s3", mass=12.12, mass_unit="g")

    # Enzyme stock
    stocksolution_dict["sto2"] = StockSolution(name="P450 stock")
    stocksolution_dict["sto2"].addComponent("s0", volume=0.001, volume_unit="l")
    stocksolution_dict["sto2"].addComponent("p0", mass=0.000008, mass_unit="g")


    all_dicts = {
        **reactant_dict,
        **enzyme_dict,
        **stocksolution_dict
    }


    #print(stocksolution._getSpecies("s1"))
    #print(stocksolution._getTotalVolume())


    #print(stocksolution_dict["sto2"].getConcentrations(all_dicts))

    measurement = Measurement(name="yolo")

    measurement.fromStockSolutions(
        stocksolution_dict=stocksolution_dict,
        all_dict=all_dicts,
        components={
            "sto0": 50,
            "sto1": 50,
            "sto2": 100
        }
    )

    print(measurement)

