#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include "atom.h"
#include "geo.h"

namespace hw1
{

    void read_geoin(hw2::geo &atom_group)
    {
        std::ifstream reader;
        std::string infoline;
        reader.open(atom_group.getInputFilePath(), std::ios::in);
        if (!reader.is_open())
        {
            std::cout << "Not find files：" << atom_group.getInputFilePath() << std::endl;
            return;
        }

        double xx, yy, zz;
        std::string useless;
        if (getline(reader, infoline) && infoline == "%CELL_PARAMETER")
        {
            reader >> xx >> useless >> useless;
            reader >> useless >> yy >> useless;
            reader >> useless >> useless >> zz;
            atom_group.setUnitsize(xx);
        }

        while (infoline != "%ATOMIC_POSTION")
        {
            getline(reader, infoline);
        }

        if (infoline == "%ATOMIC_POSTION")
        {
            hw2::atom ainfo;
            for (int i = 0; i < atom_group.get_atom_number(); i++)
            {
                reader >> useless >> xx >> yy >> zz;
                ainfo.setPostion(hw2::setXToUnit(xx, atom_group.getUnitsize()),
                                 hw2::setXToUnit(yy, atom_group.getUnitsize()),
                                 hw2::setXToUnit(zz, atom_group.getUnitsize()));
                atom_group.getAtoms().push_back(ainfo);
            }
        }

        std::cout << infoline << std::endl;

        reader.close();
        return;
    }

    void write_file(std::string filename, std::vector<hw2::atom> &atomic_information)
    {
        std::cout << atomic_information.size() << std::endl;
        std::ofstream writer;
        writer.open(filename, std::ios::out);
        if (!writer.is_open())
        {
            std::cout << "Not find files：" << filename << std::endl;
            return;
        }
        writer << "%ATOMIC_POSTION (BOHR)" << std::endl;
        for (int i = 0; i < atomic_information.size(); i++)
        {
            writer << i << "\t" << std::fixed << std::setprecision(12) << atomic_information[i].getPostionX() << "\t" << atomic_information[i].getPostionY() << "\t" << atomic_information[i].getPostionZ() << std::endl;
        }
        writer.close();
        return;
    }

} // namespace hw1

namespace hw2
{

    void read_mdin(hw2::geo &atom_group)
    {
        std::ifstream reader;
        reader.open("md.in", std::ios::in);

        std::string word;
        double info;

        reader >> word >> word;
        atom_group.setAtomname(word);
        reader >> word >> info;
        atom_group.set_atom_number((int)info);
        reader >> word >> word;
        atom_group.setInputFilePath(word);
        reader >> word >> info;
        atom_group.setRcut(info);
        reader >> word >> info;
        atom_group.setNeighbor_r(info);
        reader >> word >> info;
        atom_group.setNeighbor_n(info);
    }

    void write_No12_neighborAtom_table(std::string filename, hw2::geo &geo)
    {
        // std::cout << geo.neighborAtom_number[11] << std::endl;
        std::ofstream writer;
        writer.open(filename, std::ios::out);
        if (!writer.is_open())
        {
            std::cout << "Not find files：" << filename << std::endl;
            return;
        }
        for (int i = 0; i < geo.get_neighborAtomNumber(11); i++)
        {
            writer << geo.get_neighborAtomTable_pointer()[11][i] + 1 << "\t"
                   << std::fixed << std::setprecision(12)
                   << geo.getAtoms()[geo.get_neighborAtomTable_pointer()[11][i]].getPostionX() << "\t"
                   << geo.getAtoms()[geo.get_neighborAtomTable_pointer()[11][i]].getPostionY() << "\t"
                   << geo.getAtoms()[geo.get_neighborAtomTable_pointer()[11][i]].getPostionZ() << std::endl;
            // std::cout << geo.neighborAtom_table[11][i] << std::endl;
        }
        writer.close();
    }

} // namespace hw2

int main()
{
    std::string filenameout = "geo.out";

    hw2::geo atom_group;

    hw2::read_mdin(atom_group);

    hw1::read_geoin(atom_group);

    atom_group.init_neighborAtom_table();

    hw2::write_No12_neighborAtom_table(filenameout, atom_group);
}
