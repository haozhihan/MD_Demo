#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include "atom.h"
#include "geo.h"


int main()
{

    hw2::geo atom_group;

    atom_group.readMDIN();
    atom_group.readGeoIN();

    // std::cout << atom_group.getAtoms()[0] << std::endl;

    atom_group.init_neighborAtom_table();

    // double totalE = atom_group.total_energy();
    // atom_group.total_force();

    // hw3::write_hw3(totalE, atom_group.getAtoms());
}
