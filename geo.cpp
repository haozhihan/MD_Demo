#include "geo.h"
#include "atom.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

namespace hw2
{

    void geo::output(int step ){
        std::ofstream ofs_postion("./OUTPUT/position.txt", std::ios::app);
        std::ofstream ofs_velocity("./OUTPUT/velocity.txt", std::ios::app);
        std::ofstream ofs_force("./OUTPUT/force.txt", std::ios::app);
        std::ofstream ofs_log("./OUTPUT/run.log", std::ios::app);
        
        // postion.txt
        ofs_postion << "Step: " << step << ", Time: " << step * MDstep_time << "ps" << std::endl;
        for (int i = 0; i < atom_number; i++)
        {
            ofs_postion << i << "\t" << std::fixed << std::setprecision(12) << atoms[i].getPostionX() << "\t" << atoms[i].getPostionY() << "\t" << atoms[i].getPostionZ() <<std::endl;		
        }

        // velocity.txt
        ofs_velocity << "Step: " << step << ", Time: " << step * MDstep_time << "ps" << std::endl;
        for (int i = 0; i < atom_number; i++)
        {
            ofs_velocity << i << "\t" << std::fixed << std::setprecision(12) << atoms[i].getVelocityX() << "\t" << atoms[i].getVelocityY() << "\t" << atoms[i].getVelocityZ() <<std::endl;		
        }

        // force.txt
        ofs_force << "Step: " << step << ", Time: " << step * MDstep_time << "ps" << std::endl;
        for (int i = 0; i < atom_number; i++)
        {
            ofs_force << i << "\t" << std::fixed << std::setprecision(12) << atoms[i].getForceX() << "\t" << atoms[i].getForceY() << "\t" << atoms[i].getForceZ() <<std::endl;		
        }


        // run.log
        double total_potential = total_energy();

        double total_kinetic = 0.0;
        for (int i = 0; i < atom_number; i++)
        {
            total_kinetic = total_kinetic + pow(atoms[i].getVelocityX(), 2) + pow(atoms[i].getVelocityY(), 2) + pow(atoms[i].getVelocityZ(), 2);
        }
        total_kinetic = total_kinetic * mass / 2;

        double temperature =  2 * total_kinetic / atom_number / (1.38064852e-23 / 1.602176634e-19)   / 3 ;

        if (step == 0)
        {
            ofs_log << "Step: " << "\t" << "Time: "<< "\t" << "total_energy: "<< "\t" << "total_potential: " << "\t"  << "total_kinetic_energy: " << "\t" << "temperature: "  << std::endl;
        }
        ofs_log << step  << "\t" << std::fixed << std::setprecision(2)  << step * MDstep_time << "ps"  << "\t"  << std::fixed << std::setprecision(12) 
        << total_potential + total_kinetic << "\t"  << total_potential  << "\t" << total_kinetic  << "\t" << temperature << std::endl;
    }

    void geo::runMD()
    {

        // define store pre force value
        double** pre_atoms_force_matrix = new double* [atom_number];
        for (int i = 0; i < atom_number; i++)
        {
            pre_atoms_force_matrix[i] = new double[3];
            pre_atoms_force_matrix[i][0] = 0;
            pre_atoms_force_matrix[i][1] = 0;
            pre_atoms_force_matrix[i][2] = 0;
        }
        
        for (int step = 0; step <= total_step; step++)
        {
            // output
            if (step % output_everystep == 0)
            {
                output(step);
            }

            // update neighborAtom_table
            if (step % update_neighbor_step == 0 && step != 0)
            {
                build_neighborAtom_table();
            }

            // store pre force
            for (int i = 0; i < atom_number; i++)
            {
                pre_atoms_force_matrix[i][0] = atoms[i].getForceX();
                pre_atoms_force_matrix[i][1] = atoms[i].getForceY();
                pre_atoms_force_matrix[i][2] = atoms[i].getForceZ();
            }

            // update postion
            for (int i = 0; i < atom_number; i++)
            {
                double x = atoms[i].getPostionX() + atoms[i].getVelocityX() * MDstep_time + atoms[i].getForceX() * pow(MDstep_time, 2) / 2 / mass;
                double y = atoms[i].getPostionY() + atoms[i].getVelocityY() * MDstep_time + atoms[i].getForceY() * pow(MDstep_time, 2) / 2 / mass;
                double z = atoms[i].getPostionZ() + atoms[i].getVelocityZ() * MDstep_time + atoms[i].getForceZ() * pow(MDstep_time, 2) / 2 / mass;
                atoms[i].setPostion(x, y, z);           
            }

            // cal new force according new postion
            cal_every_atom_force();

            // update velocity
            for (int i = 0; i < atom_number; i++)
            {
                double x = atoms[i].getVelocityX() + MDstep_time * ( atoms[i].getForceX() + pre_atoms_force_matrix[i][0] ) / 2 / mass;
                double y = atoms[i].getVelocityY() + MDstep_time * ( atoms[i].getForceY() + pre_atoms_force_matrix[i][1] ) / 2 / mass;
                double z = atoms[i].getVelocityZ() + MDstep_time * ( atoms[i].getForceZ() + pre_atoms_force_matrix[i][2] ) / 2 / mass;
                atoms[i].setVelocity(x, y, z);
            }

            
        }
        
    }

    void geo::readMDIN(std::string mdINpath)
    {
        std::ifstream reader;
        reader.open(mdINpath, std::ios::in);
        
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


        reader >> word >> word;
        this->is_read_velocity = word;

        // important!!!
        reader >> word >> info;
        this->mass = info * 10 / 6.02214086e23 / 1.602176634e-19;

        reader >> word >> info;
        this->update_neighbor_step = (int)info;

        reader >> word >> info;
        this->total_step = (int)info;

        reader >> word >> info;
        this->MDstep_time = info;

        reader >> word >> info;
        this->output_everystep = info;
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


        if (this->is_read_velocity == "yes")
        {
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


    void geo::build_neighborAtom_table()
    {

        // init
        for (int i = 0; i < atom_number; i++)
        {
            neighborAtom_number[i] = 0;
            for (int j = 0; j < neighbor_n; j++)
            {
                neighborAtom_table[i][j] = 0;
            }
        }

        // calculate
        for (int i = 0; i < atom_number; i++)
        {
            for (int j = 0; j < atom_number; j++)
            {
                if (i == j)
                {
                    continue;
                }
                double dist = calculate_distance(atoms[i].getPostionX(), atoms[i].getPostionY(), atoms[i].getPostionZ(),
                                                 atoms[j].getPostionX(), atoms[j].getPostionY(), atoms[j].getPostionZ(), unit_size);

                if (dist < get_R_neighborAtom())
                {
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

    void geo::cal_every_atom_force()
    {   
        for (int i = 0; i < atom_number; i++)
        {
            atoms[i].setForceZero();
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
