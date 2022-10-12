#include "geo.h"
#include "atom.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

namespace hw2
{

    void geo::readMDIN()
    {
        std::ifstream reader;
        reader.open("md.in", std::ios::in);
        
        std::string word;
        double info;

        reader >> word >> word;
        this->atom_name = word;

        reader >> word >> info;
        this->atom_number = (int)info;

        reader >> word >> word;
        this->geoin_file_path = word;

        reader >> word >> info;
        this->rcut = info;

        reader >> word >> info;
        this->neighbor_r = info;

        reader >> word >> info;
        this->neighbor_n = info;

        reader >> word >> info;
        this->epsilon = info;

        reader >> word >> info;
        this->sigma = info;

        this->ecut = 4 * this->epsilon * (pow(this->sigma / rcut, 12) - pow(this->sigma / rcut, 6));
    }



    void geo::readGeoIN()
    {
        std::ifstream reader;
        std::string infoline;
        reader.open(geoin_file_path, std::ios::in);
        if (!reader.is_open())
        {
            std::cout << "Not find files：" << geoin_file_path << std::endl;
            return;
        }
        double xx, yy, zz;
        std::string useless;

        // read %CELL_PARAMETER
        if (getline(reader, infoline) && infoline == "%CELL_PARAMETER")
        {
            reader >> xx >> useless >> useless;
            reader >> useless >> yy >> useless;
            reader >> useless >> useless >> zz;
            this->unit_size = xx;
        }


        // read ATOMIC_POSTION
        while (infoline != "%ATOMIC_POSTION")
        {
            getline(reader, infoline);
        }

        if (infoline == "%ATOMIC_POSTION")
        {
            hw2::atom ainfo;
            for (int i = 0; i < atom_number; i++)
            {
                reader >> useless >> xx >> yy >> zz;
                ainfo.setPostion(hw2::setXToUnit(xx, this->unit_size),
                                 hw2::setXToUnit(yy, this->unit_size),
                                 hw2::setXToUnit(zz, this->unit_size));
                atoms.push_back(ainfo);
            }
        }

        // read ATOMIC_VELOCITY
        while (infoline != "%ATOMIC_VELOCITY")
        {
            getline(reader, infoline);
        }

        if (infoline == "%ATOMIC_VELOCITY")
        {
            for (int i = 0; i < atom_number; i++)
            {
                reader >> useless >> xx >> yy >> zz;
                atoms[i].setVelocity(xx, yy, zz);
            }
        }

        reader.close();
        return;
    }


    std::vector<atom> &geo::getAtoms()
    {
        return this->atoms;
    }

    double geo::get_R_neighborAtom()
    {
        return rcut + neighbor_r;
    }


    void geo::init_neighborAtom_table()
    {

        // init
        neighborAtom_number = new int[atom_number]();
        neighborAtom_table = new int *[atom_number]();
        for (int i = 0; i < atom_number; i++)
        {
            neighborAtom_table[i] = new int[neighbor_n];
            neighborAtom_number[i] = 0;
            for (int j = 0; j < neighbor_n; j++)
            {
                neighborAtom_table[i][j] = 0;
            }
        }

        // std::cout << atomic_information.size() << std::endl;

        // calculate
        for (int i = 0; i < atom_number; i++)
        {

            // std::cout << "/* message */" << std::endl;
            for (int j = 0; j < atom_number; j++)
            {
                if (i == j)
                {
                    continue;
                }
                double dist = calculate_distance(atoms[i].getPostionX(), atoms[i].getPostionY(), atoms[i].getPostionZ(),
                                                 atoms[j].getPostionX(), atoms[j].getPostionY(), atoms[j].getPostionZ(), unit_size);

                // std::cout << distance << std::endl;

                if (dist < get_R_neighborAtom())
                {
                    // if(neighborAtom_number[i] >= neighbor_n)
                    // {

                    // }

                    // std::cout << "/* message */" << std::endl;

                    if (neighborAtom_number[i] < neighbor_n)
                    {
                        neighborAtom_table[i][neighborAtom_number[i]] = j;
                        neighborAtom_number[i]++;
                    }
                }
            }
        }
    }


    int geo::get_neighborAtomNumber(int n)
    {
        return this->neighborAtom_number[n];
    }

    int **geo::get_neighborAtomTable_pointer()
    {
        return this->neighborAtom_table;
    }

    double geo::total_energy()
    {
        double E = 0.0;
        for (int i = 0; i < atom_number; i++)
        {
            for (int j = 0; j < neighborAtom_number[i]; j++)
            {
                int atomnumber = neighborAtom_table[i][j];
                E = E + atoms[i].twoAtompPotential(atoms[atomnumber], unit_size, rcut, epsilon, sigma, ecut);
            }
        }
        return E / 2;
    }

    void geo::total_force()
    {
        for (int i = 0; i < atom_number; i++)
        {
            for (int j = 0; j < neighborAtom_number[i]; j++)
            {
                int atomnumber = neighborAtom_table[i][j];
                atoms[i].twoAtomForce(atoms[atomnumber], unit_size, rcut, epsilon, sigma);
            }
        }
    }

    double calculate_distance(double x1, double y1, double z1, double x2, double y2, double z2, double size)
    {
        double distance_x = short_distance_OneDirection(x1, x2, size);
        double distance_y = short_distance_OneDirection(y1, y2, size);
        double distance_z = short_distance_OneDirection(z1, z2, size);
        return std::sqrt(std::pow(distance_x, 2) + std::pow(distance_y, 2) + std::pow(distance_z, 2));
    }

    double short_distance_OneDirection(double x, double y, double size)
    {
        // double one = abs(x - y);
        // double two = abs(x - y - size);
        // double three = abs(x + size - y);

        if (x < y)
        {
            double s1 = y - x;
            double s2 = size - y + x;
            if (s1 < s2)
            {
                return x - y;
            }
            else
            {
                return s2;
            }
        }
        else
        {
            double s1 = x - y;
            double s2 = size - x + y;
            if (s1 < s2)
            {
                return s1;
            }
            else
            {
                return x - y - size;
            }
        }

        // return std::min(std::min(one, two), three);
    }

} // namespace hw2
