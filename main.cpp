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

    // std::cout << atom_group.total_energy() << std::endl;


    atom_group.runMD();
    // double totalE = atom_group.total_energy();
    // atom_group.cal_every_atom_force();
}
