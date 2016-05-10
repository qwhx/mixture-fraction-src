/* 
This is a small library class to calculate mixture fraction from OpenFOAM
mixture class. 

Author: Jian Cai 
Contact: jcai@uwyo.edu
Please submit comments/questions/requests directly to the source code 
host website.
*/

#include "mixtureFraction.H"

Foam::mixtureFraction::mixtureFraction
(
    const Foam::basicMultiComponentMixture& comp,
    const Foam::fvMesh& mesh
)
:
composition_(comp),
mesh_(mesh),
definition_
(
    IOobject
    (
        "mixtureFractionProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
speciesData_
(
    IOobject
    (
        "elementSpeciesProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
mf_
(
    IOobject
    (
        definition_.lookup("mixtureFractionName"),
        mesh_.time().timeName(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("mixFrac", dimless, 0.0)
),
elements_(),
elementWeights_(),
species_(),
fuel_(),
oxy_(),
fuelNorm_(0.0),
oxyNorm_(0.0)
{
    parseDict();
    calculate();
}

Foam::mixtureFraction::~mixtureFraction()
{}

const Foam::volScalarField& Foam::mixtureFraction::mixFrac()
{
    return mf_;
}

void Foam::mixtureFraction::calculate()
{
    mf_=0.0;
    forAll (species_, isp)
    {
        const volScalarField& Y=composition_.Y(isp);
        const scalarField& Ye = species_[isp];
        forAll (elementWeights_, ie)
        {
            const scalar ew = elementWeights_[ie];
            if (ew<=0)
            {
                continue;
            }

            mf_ += Y * Ye[ie] * ew;
        }
    }
    mf_ = (mf_-oxyNorm_)/(fuelNorm_-oxyNorm_);
}

void Foam::mixtureFraction::parseDict()
{
    // parse elements
    wordList elemNames (speciesData_.lookup("elements"));
    elements_ = speciesTable(elemNames);

    // parse mixture fraction definition
    elementWeights_.resize(elements_.size(), 0.0);
    const dictionary & weightDict = definition_.subDict("elementWeights");
    forAllConstIter(IDLList<entry>, weightDict, iter)
    {
        const word elementName = iter().keyword();
        if (!elements_.contains(elementName))
        {
            FatalErrorIn
            (
                "mixtureFraction::parseDict"
            )   << "Weight factor of element "<< elementName
                << " is specified, but this element is not contained "
                << " in elementSpeciesProperties dictionary. " << nl
                << exit(FatalError);
        }
        elementWeights_[elements_[elementName]] = 
                    pTraits<scalar>(iter().stream());
    }

    // parse species information
    const speciesTable sp (composition_.species());
    const label nsp = sp.size();

    species_.setSize(nsp);

    forAll (sp, isp)
    {
        const word speciesName = sp[isp];
        if (
                (!speciesData_.found(speciesName))
              ||(!speciesData_.isDict(speciesName))
            )
        {
            FatalErrorIn
            (
                "mixtureFraction::parseDict"
            )   << "Property subdictionary is not found for species "
                << speciesName << nl << exit(FatalError);
        }

        const dictionary & spData = speciesData_.subDict(speciesName);
        species_[isp].resize(elements_.size(), 0.0);
        forAll (elements_, ie)
        {
            species_[isp][ie] = spData.template lookupOrDefault<scalar>
                (elements_[ie], 0.0);
        }
    }

    // parse fuel stream properties
    const dictionary & fuel = definition_.subDict("fuelStream");
    fuel_.resize(sp.size(), 0.0);
    fuelNorm_=0.0;
    forAll (sp, isp)
    {
        fuel_[isp] = fuel.template lookupOrDefault<scalar>(sp[isp], 0.0);
        fuelNorm_+=sum(elementWeights_*species_[isp])*fuel_[isp];
    }

    // parse oxydizer stream properties
    const dictionary & oxy = definition_.subDict("oxydizerStream");
    oxy_.resize(sp.size(), 0.0);
    oxyNorm_ = 0.0;
    forAll (sp, isp)
    {
        oxy_[isp] = oxy.template lookupOrDefault<scalar>(sp[isp], 0.0);
        oxyNorm_+=sum(elementWeights_*species_[isp])*oxy_[isp];
    }
}

