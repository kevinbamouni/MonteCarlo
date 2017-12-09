#include <random>
#include <iostream>

double alea_u();

int main()
{
    double a(alea_u());
    std::cout << a;
}

double alea_u(){
    // generate random uniform numbers by c++11
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    //Use dis to transform the random unsigned int generated by gen into a double in [1, 2)
    return dis(gen); //Each call to dis(gen) generates a new random double
}
