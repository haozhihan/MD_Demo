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

    hw2::geo atom_group("./INPUT/md.in");

    // hw2::randomNumber(-0.5, 0.5);

    std::cout << hw2::randomNumber(-0.5, 0.5) << std::endl;
    std::cout << hw2::randomNumber(-0.5, 0.5) << std::endl;
    std::cout << hw2::randomNumber(-0.5, 0.5) << std::endl;

    // std::cout << atom_group.getAtoms()[0].getVelocityX() << std::endl; 
    
    // atom_group.runMD();

}
