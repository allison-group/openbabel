/**********************************************************************
Copyright (C) 2002-2006 by Geoffrey Hutchison
Some portions Copyright (C) 2004 by Chris Morley

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include <openbabel/babelconfig.h>

#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>


using namespace std;
namespace OpenBabel
{

  class PSI4Format : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PSI4Format()
    {
      OBConversion::RegisterFormat("psi",this);
      OBConversion::RegisterFormat("psi4",this);
    }

    virtual const char* Description() //required
    {
      return
        "PSI4 input format\n"
        "The input format for the quantum-mechanics program PSI4.\n"
        "Write Options e.g. -xc\n"
        "  c  Write an input file for the CNDO/INDO program.\n\n";
    };

    virtual const char* SpecificationURL()
    {return "";}; //optional

    //Flags() can return be any the following combined by | or be omitted if none apply
    // NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY
    virtual unsigned int Flags()
    {
      return NOTREADABLE | WRITEONEONLY;
    };

    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };

  //Make an instance of the format class
  PSI4Format thePSI4Format;

  ////////////////////////////////////////////////////////////////

  bool PSI4Format::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if (pmol == nullptr)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    ofs << "molecule " << pmol->GetTitle() << " {\n";
    ofs << "  " << pmol->GetTotalCharge() << "  " << pmol->GetTotalSpinMultiplicity() << '\n';
    ofs << "# Unoptimised Energy: " << pmol->GetEnergy() << '\n';

    char buffer[BUFF_SIZE];
    vector<OBAtom*>::iterator i;
    OBAtom *atom;

    for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i))
      {
        snprintf(buffer, BUFF_SIZE, "  %2s   %10.6f %10.6f %10.6f",
                 OBElements::GetSymbol(atom->GetAtomicNum()),
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ());
        ofs << buffer << '\n';
      }
    ofs << "  units angstrom\n";
    ofs << "}\n";

    bool nonRestricted = pmol->GetTotalSpinMultiplicity() == 1;
    ofs << "set basis 6-31G(d)\n"
        << "set reference ";
    if (nonRestricted) { ofs << "rhf\n"; } else { ofs << "uhf\n"; }
    ofs << "set opt_coordinates both\n"
        << "energy, wfn = optimize(\"b3lyp\", return_wfn = True)\n"
        << "fchk(wfn, \"" << pmol->GetTitle() << ".fchk\")\n";



    return(true);
  }

} //namespace OpenBabel
