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
    
    atom_group.runMD();
}
