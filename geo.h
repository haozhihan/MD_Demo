#include <string>
#include <vector>
#include "atom.h"

namespace hw2
{
    class geo
    {
    private:
        // from geo1.in
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

        // hw4 add
        std::string is_read_velocity;
        double mass;
        int update_neighbor_step;
        int total_step;
        double MDstep_time;
        int output_everystep;

    public:
        geo(std::string mdINpath){
            readMDIN(mdINpath);
            readGeoIN();
                    
            // init 0 for neighborAtom_number & neighborAtom_table 
            neighborAtom_number = new int[atom_number]();
            neighborAtom_table = new int *[atom_number]();
            for (int i = 0; i < atom_number; i++)
            {
                neighborAtom_table[i] = new int[neighbor_n];
            }

            build_neighborAtom_table();


            // cal init force for every atom
            cal_every_atom_force();
        };

        ~geo(){};


        void readMDIN(std::string mdINpath);
        void readGeoIN();

        std::vector<atom> &getAtoms();
        double get_R_neighborAtom();
        

        void build_neighborAtom_table();

        int get_neighborAtomNumber(int n);

        int **get_neighborAtomTable_pointer();

        double total_energy();

        void cal_every_atom_force();

        void runMD();

        void output(int step);
    };

    double calculate_distance(double x1, double y1, double z1, double x2, double y2, double z2, double size);

    double short_distance_OneDirection(double x, double y, double size);

} // namespace hw2
