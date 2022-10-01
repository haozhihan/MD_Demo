#include <string>
#include <vector>
#include "atom.h"

namespace hw2
{
    class geo
    {
    private:
        double sizeX;
        double sizeY;
        double sizeZ;

        int atom_number;
        std::string atom_type;



    public:
        geo(/* args */){};
        ~geo(){};

        double rcut;
        double neighbor_r;
        int neighbor_n;

        int* neighborAtom_number;
        int** neighborAtom_table;

        void set_sizeX(double x);
        void set_sizeY(double y);
        void set_sizeZ(double z);

        double get_sizeX();
        double get_sizeY();
        double get_sizeZ();


        void set_atom_number(int number);

        int get_atom_number();

        double get_neighbor_rcut();

        // void init_neighborAtom_number(std::vector<hw2::atom>& atomic_information);
        
        void init_neighborAtom_table(std::vector<hw2::atom>& atomic_information,  double distance, double size);

        // void get_neighborAtom_table(){
        //     // return neighborAtom_table
        // }
    };

    double calculate_distance(atom atom1info, atom atom2info, double size);

    double short_distance_OneDirection(double x, double y, double size);

    // geo::geo(/* args */)
    // {
    // }

    // geo::~geo()
    // {
    // }
    
} // namespace hw2


