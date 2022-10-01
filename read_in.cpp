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

    void read_file(std::string filename, std::vector<hw2::atom>& atomic_information, hw2::geo &atom_group){
        std::ifstream reader;
        std::string info;
        reader.open(filename, std::ios::in);
        if(!reader.is_open()){
            std::cout << "Not find files：" << filename << std::endl;
            return;
        }


        double xx,yy,zz,useless;        
        if (getline(reader, info) && info == "%CELL_PARAMETER"){
            reader >> xx >> useless >> useless;
            reader >> useless >> yy >> useless;
            reader >> useless >> useless >> zz;
            atom_group.set_sizeX(xx);
            atom_group.set_sizeY(yy);
            atom_group.set_sizeZ(zz);
        }
        
        while (info != "%ATOMIC_POSTION"){
            getline(reader, info);
        }
        std::string he,x,y,z;
        while( getline(reader, info) && !info.empty()){
            std::stringstream word(info);
            word >> he;
            word >> x;
            word >> y;
            word >> z;
            hw2::atom ainfo;
            ainfo.setPostion( hw2::setXToUnit(atof(x.c_str()), atom_group.get_sizeX()), 
                              hw2::setXToUnit(atof(y.c_str()), atom_group.get_sizeY()), 
                              hw2::setXToUnit(atof(z.c_str()), atom_group.get_sizeZ()) );

            // ainfo.setPostion( atof(x.c_str()), 
            //                   atof(y.c_str()), 
            //                   atof(z.c_str()) );
            atomic_information.push_back(ainfo);
        }
        reader.close();
        return;
    }

    void write_file(std::string filename, std::vector<hw2::atom>& atomic_information){
        std::cout << atomic_information.size() << std::endl;
        std::ofstream writer;
        writer.open(filename, std::ios::out);
        if(!writer.is_open()){
            std::cout << "Not find files：" << filename << std::endl;
            return;
        }
        writer << "%ATOMIC_POSTION (BOHR)" << std::endl;
        for (int i = 0; i < atomic_information.size(); i++){
            writer << i << "\t" << std::fixed << std::setprecision(12) << atomic_information[i].getPostionX() << "\t" << atomic_information[i].getPostionY() << "\t" << atomic_information[i].getPostionZ() << std::endl;
        }
        writer.close();
        return;
    }
    
} // namespace hw1


namespace hw2
{
    void write_No12_neighborAtom_table(std::string filename, hw2::geo &geo, std::vector<hw2::atom>& atomic_information){
        // std::cout << geo.neighborAtom_number[11] << std::endl;
        std::ofstream writer;
        writer.open(filename, std::ios::out);
        if(!writer.is_open()){
            std::cout << "Not find files：" << filename << std::endl;
            return;
        }
        for (int i = 0; i < geo.neighborAtom_number[11]; i++)
        {
            writer << geo.neighborAtom_table[11][i] + 1 << "\t" << std::fixed << std::setprecision(12) << atomic_information[geo.neighborAtom_table[11][i]].getPostionX() << "\t" << atomic_information[geo.neighborAtom_table[11][i]].getPostionY() << "\t" << atomic_information[geo.neighborAtom_table[11][i]].getPostionZ() << std::endl;
            // std::cout << geo.neighborAtom_table[11][i] << std::endl;
        }
        writer.close();
    }
} // namespace hw2



int main(){
    std::string filenamein;
    std::string filenameout = "geo.out";
    std::vector<hw2::atom> atomic_information;
    hw2::geo atom_group;

    std::string useless;
    double info;
    std::ifstream reader;
    reader.open("md.in", std::ios::in);
    reader >> useless >> useless;
    reader >> useless >> useless;
    reader >> useless >> filenamein;
    reader >> useless >> info;
    atom_group.rcut = info;
    reader >> useless >> info;
    atom_group.neighbor_r = info;
    reader >> useless >> info;
    atom_group.neighbor_n = info;
    


    hw1::read_file(filenamein, atomic_information, atom_group);

    // std::cout << atomic_information.size() << std::endl;

    // hw1::write_file(filenameout, atomic_information);

    atom_group.set_atom_number(atomic_information.size());

    atom_group.init_neighborAtom_table(atomic_information, atom_group.get_neighbor_rcut() ,atom_group.get_sizeX());

    hw2::write_No12_neighborAtom_table(filenameout, atom_group, atomic_information);

}



