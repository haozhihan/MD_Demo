#include "geo.h"
#include "atom.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace hw2
{

    double geo::get_neighbor_rcut(){
        return rcut + neighbor_r;
    }

    void geo::set_sizeX(double x)
    {
        sizeX = x;
    }

    void geo::set_sizeY(double y)
    {
        sizeY = y;
    }

    void geo::set_sizeZ(double z)
    {
        sizeZ = z;
    }

    double geo::get_sizeX()
    {
        return sizeX;
    }

    double geo::get_sizeY()
    {
        return sizeY;
    }

    double geo::get_sizeZ()
    {
        return sizeZ;
    }

    int geo::get_atom_number()
    {
        return atom_number;
    }

    void geo::set_atom_number(int number)
    {
        atom_number = number;
    }

    void geo::init_neighborAtom_table(std::vector<hw2::atom> &atomic_information, double distance, double size)
    {

        // init
        neighborAtom_number = new int[atomic_information.size()]();
        neighborAtom_table = new int *[atomic_information.size()]();
        for (int i = 0; i < atomic_information.size(); i++)
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
        for (int i = 0; i < atomic_information.size(); i++)
        {

            // std::cout << "/* message */" << std::endl;
            for (int j = 0; j < atomic_information.size(); j++)
            {
                if (i == j)
                {
                    continue;
                }
                double dist = calculate_distance(atomic_information[i], atomic_information[j], size);

                // std::cout << distance << std::endl;
                
                if (dist < distance)
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

    double calculate_distance(atom atom1info, atom atom2info, double size)
    {
        double distance_x = short_distance_OneDirection(atom1info.getPostionX(), atom2info.getPostionX(), size);
        double distance_y = short_distance_OneDirection(atom1info.getPostionY(), atom2info.getPostionY(), size);
        double distance_z = short_distance_OneDirection(atom1info.getPostionZ(), atom2info.getPostionZ(), size);
        return std::sqrt(std::pow(distance_x, 2) + std::pow(distance_y, 2) + std::pow(distance_z, 2));
    }

    double short_distance_OneDirection(double x, double y, double size)
    {
        // double one = abs(x - y);
        // double two = abs(x - y - size);
        // double three = abs(x + size - y);

        if (x < y)
        {
            double s1 = y -x;
            double s2 = size - y + x;
            if (s1 < s2)
            {
                return s1;
            } else
            {
                return s2;
            }
        } else {
            double s1 = x - y;
            double s2 = size - x + y;
            if (s1 < s2)
            {
                return s1;
            } else
            {
                return s2;
            }
        }
        

        // return std::min(std::min(one, two), three);
    }

} // namespace hw2
