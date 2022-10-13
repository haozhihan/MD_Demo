#ifndef ATOM_H
#define ATOM_H

#include <string>

namespace hw2
{

    class atom
    {
    private:
        std::string type;
        double postion_x;
        double postion_y;
        double postion_z;

        double velocity_x;
        double velocity_y;
        double velocity_z;

        double force_x;
        double force_y;
        double force_z;

    public:
        atom(/* args */){};
        ~atom(){};

        void setPostion(double in_postion_x, double in_postion_y, double in_postion_z);

        void setVelocity(double in_velocity_x, double in_velocity_y, double in_velocity_z);


        void setForceZero();

        double getPostionX();

        double getPostionY();

        double getPostionZ();


        double getVelocityX();

        double getVelocityY();

        double getVelocityZ();


        double getForceX();

        double getForceY();

        double getForceZ();

        double twoAtompPotential(atom another, double size, double rcut, double epsilon, double sigma, double ecut);

        void twoAtomForce(atom another, double size, double rcut, double epsilon, double sigma);
    };

    double setXToUnit(double x, double sizeX);

    // atom::atom(/* args */)
    // {
    // }

    // atom::~atom()
    // {
    // }

}

#endif