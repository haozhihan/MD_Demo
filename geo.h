#include <string>
#include <vector>
#include "atom.h"

namespace hw2
{
    class geo
    {
    private:
        double unit_size;
        std::vector<atom> atoms;

        int *neighborAtom_number;
        int **neighborAtom_table;

        // md.in
        std::string atom_name;
        int atom_number;
        std::string geoin_file_path;

        double rcut;
        double neighbor_r;
        int neighbor_n;

// hw3 add
        double ecut;
        double epsilon;
        double sigma;

    public:
        geo(/* args */){};
        ~geo(){};


        void readMDIN();
        void readGeoIN();

        std::vector<atom> &getAtoms();
        double get_R_neighborAtom();
        

        void init_neighborAtom_table();

        int get_neighborAtomNumber(int n);

        int **get_neighborAtomTable_pointer();

        double total_energy();

        void total_force();
    };

    double calculate_distance(double x1, double y1, double z1, double x2, double y2, double z2, double size);

    double short_distance_OneDirection(double x, double y, double size);

} // namespace hw2
