#include "atom.h"
#include "geo.h"
#include <cmath>
#include <random>

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

    double randomNumber(double x1, double x2)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(x1, x2);
        return dis(gen);
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

    void atom::setForceZero(){
        this->force_x = 0;
        this->force_y = 0;
        this->force_z = 0;
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

    double atom::getVelocityX()
    {
        return this->velocity_x;
    }
    double atom::getVelocityY()
    {
        return this->velocity_y;
    }
    double atom::getVelocityZ()
    {
        return this->velocity_z;
    }


    double atom::getForceX(){
        return this->force_x;
    }

    double atom::getForceY(){
        return this->force_y;
    }

    double atom::getForceZ(){
        return this->force_z;
    }

    double atom::twoAtompPotential(atom another, double size, double rcut, double epsilon, double sigma, double ecut){
        
        double r = hw2::calculate_distance( this->getPostionX(), this->getPostionY(), this->getPostionZ(), 
                                            another.getPostionX(), another.getPostionY(), another.getPostionZ(), size);
        if (r <= rcut)
        {
            double E = 4 * epsilon * (pow(sigma/r, 12) - pow(sigma/r, 6)) - ecut;
            return E;
        } else 
        {
            return 0.0;
        }

    }

    void atom::twoAtomForce(atom another, double size, double rcut, double epsilon, double sigma){

        double r = hw2::calculate_distance( this->getPostionX(), this->getPostionY(), this->getPostionZ(), 
                                            another.getPostionX(), another.getPostionY(), another.getPostionZ(), size);
        if (r <= rcut)
        {
            double temp = 4 * epsilon * ( 12 * pow(sigma/r, 12) - 6 * pow(sigma/r, 6) ) / pow(r, 2) ;
            force_x = force_x + temp * short_distance_OneDirection(this->getPostionX(), another.getPostionX(), size);
            force_y = force_y + temp * short_distance_OneDirection(this->getPostionY(), another.getPostionY(), size);
            force_z = force_z + temp * short_distance_OneDirection(this->getPostionZ(), another.getPostionZ(), size);
        } else
        {
            force_x = force_x + 0;
            force_y = force_y + 0;
            force_z = force_z + 0;
        }
        return;
    }

} // namespace hw2
