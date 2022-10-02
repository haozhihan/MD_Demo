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

    public:
        atom(/* args */){};
        ~atom(){};

        void setPostion(double in_postion_x, double in_postion_y, double in_postion_z);

        void setVelocity(double in_velocity_x, double in_velocity_y, double in_velocity_z);

        double getPostionX();

        double getPostionY();

        double getPostionZ();
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