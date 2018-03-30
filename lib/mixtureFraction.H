/* 
This is a small library class to calculate mixture fraction from OpenFOAM
mixture class. 
*/

#ifndef mixtureFraction_H
#define mixtureFraction_H

#include "basicMultiComponentMixture.H"

namespace Foam
{

class mixtureFraction
{

    const basicMultiComponentMixture& composition_;

    const fvMesh& mesh_;

    const IOdictionary definition_;

    const IOdictionary speciesData_;

    volScalarField mf_; // mixture fraction field

    speciesTable elements_;

    scalarField elementWeights_;

    List<scalarField> species_;

    scalarList fuel_;

    scalarList oxy_;

    scalar fuelNorm_;

    scalar oxyNorm_;

public:
    // Constructor
    mixtureFraction
    (
        const basicMultiComponentMixture&,
        const fvMesh&
    );

    // Destructor
    ~mixtureFraction();

    // Access
    const volScalarField& mixFrac();
    
    // Update
    void calculate();
private:
    // Internal functions
    void parseDict();

};

} // End namespace Foam
#endif
