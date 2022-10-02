#include "atom.h"

namespace hw2
{
    double setXToUnit(double x, double sizeX)
    {
        while (x >= sizeX)
        {
            x = x - sizeX;
        }
        while (x < 0)
        {
            x = x + sizeX;
        }
        return x;
    }

    void atom::setPostion(double in_postion_x, double in_postion_y, double in_postion_z)
    {
        postion_x = in_postion_x;
        postion_y = in_postion_y;
        postion_z = in_postion_z;
    }

    void atom::setVelocity(double in_velocity_x, double in_velocity_y, double in_velocity_z)
    {
        velocity_x = in_velocity_x;
        velocity_y = in_velocity_y;
        velocity_z = in_velocity_z;
    }

    double atom::getPostionX()
    {
        return this->postion_x;
    }

    double atom::getPostionY()
    {
        return this->postion_y;
    }

    double atom::getPostionZ()
    {
        return this->postion_z;
    }
} // namespace hw2
