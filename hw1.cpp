#include <iostream>
#include <fstream> 
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

namespace hw1
{
    struct atomic_info
    {
        // std::string name;
        double postion_x;
        double postion_y;
        double postion_z;
    };

    void read_file(std::string filename, std::vector<atomic_info>& atomic_information){
        std::ifstream reader;
        std::string info;
        reader.open(filename, std::ios::in);
        if(!reader.is_open()){
            std::cout << "无法找到文件：" << filename << std::endl;
            return;
        }
        getline(reader, info);
        while (info != "%ATOMIC_POSTION"){
            getline(reader, info);
            
        }
        std::string he,x,y,z;
        double bohr = 1.8897161646321;
        while( getline(reader, info) ){
            std::stringstream word(info);
            word >> he;
            word >> x;
            word >> y;
            word >> z;
            atomic_info ainfo;
            ainfo.postion_x = atof(x.c_str()) * bohr;
            ainfo.postion_y = atof(y.c_str()) * bohr;
            ainfo.postion_z = atof(z.c_str()) * bohr;
            atomic_information.push_back(ainfo);
        }
        reader.close();
        return;
    }

    void write_file(std::string filename, std::vector<atomic_info>& atomic_information){
        std::ofstream writer;
        writer.open(filename, std::ios::out);
        if(!writer.is_open()){
            std::cout << "无法找到文件：" << filename << std::endl;
            return;
        }
        writer << "%ATOMIC_POSTION (BOHR)" << std::endl;
        for (int i = 0; i < atomic_information.size(); i++){
            writer << i << "\t" << std::fixed << std::setprecision(12) << atomic_information[i].postion_x << "\t" << atomic_information[i].postion_y << "\t" << atomic_information[i].postion_z << std::endl;
        }
        writer.close();
        return;
    }
    
} // namespace hw1


int main(){
    std::string filenamein = "geo.in";
    std::string filenameout = "geo.out";
    std::vector<hw1::atomic_info> atomic_information;

    hw1::read_file(filenamein, atomic_information);
    hw1::write_file(filenameout, atomic_information);
}



