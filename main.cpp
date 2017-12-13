#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>


// fonctions declarations
double alea_uniform_dist();
double alea_normal_dist(const double& stddev);

// fonction de calcul de moyenne
template<typename T, typename S>
T moyenne(S tab[], int& taille){
    T somme(0);
    for (int i(0);i<taille;i++){
        somme = somme + tab[i];
    }
    return somme/taille;
}

// fonction de calcul de variance sans biais
template<typename T, typename S>
T variance( S tab[],  int& taille){
    T variance_quad(0);
    for (int i(0);i<taille;i++){
        variance_quad = variance_quad + (tab[i]-moyenne<T,S>(tab, taille))*(tab[i]-moyenne<T,S>(tab, taille));
    }
    return variance_quad/(taille-1);
}

//*********************************************************************************************************
//************************************************MAIN*****************************************************
int main()
{

    // Print Application title
    std::cout << "EUROPEAN CALL OPTION PRICING WITH THE MONTE CARLO METHOD"<<std::endl;
    std::cout << "********************************************************"<<std::endl;

    // 1 PARAMETERS INITIALISATION
    // Parameters of the call option that we want to price
    double S0(100); // actual price of the underlying
    double K(100); // the strike of the call option
    double r(0.05); // risk free rate
    double T(1.0); // the maturity of the option
    double sigma(0.10); // the volatility of the underlying, for simplification we take it constant
    int N(252); // number of period per year
    int M(10000); // number of simulation;1000-1,18sec-4.3802 10000-11,649sec-4.37739; 100000-114,556
    double S[N+1]; S[0] = S0;// the path of the underlying, the fist element is the actual price of the underlying
    double dt(T/N); // time step of the underlying simulation
    double payoffs[M]; //
    double premium(0);
    double stddev(1.0);

    // print Input parameters
    std::cout << "THE CALL PARAMETERS :"<<std::endl;
    std::cout << "S0 = " << S0 <<std::endl;
    std::cout << "K = " << K <<std::endl;
    std::cout << "r = " << r <<std::endl;
    std::cout << "T = " << T <<std::endl;
    std::cout << "sigma = " << sigma <<std::endl;
    std::cout << "Monte carlo number of simulations = " << M <<std::endl;
    std::cout << "********************************************************"<<std::endl;

    // 2 Payoffs simulations loop
    for(int j(0);j<M;j++){

        // 3 path simulation loop
        for(int i(0);i<N+1;i++){
        S[i+1] = S[i]*exp((r-(0.5*sigma*sigma))*dt+sigma*sqrt(dt)*alea_normal_dist(stddev));
        }

        // 4 compute the sum of payoffs of the simulations
        payoffs[j] = std::max(S[N] - K,0.0);
    }

    // 5 Dicounted expected premium computation
    double moypayoff(moyenne<double,double>(payoffs,M)); // simulated payoff mean estimation
    premium = exp(-r*T)*(moypayoff); // Actualise la moyenne des payoff pour fournir le prix

    // estimation precisions details
    double pay_stddev(0); // standard deviation of the payoffs simulated
    pay_stddev = sqrt(variance<double,double>(payoffs,M)); // standard deviation of the payoffs simulated

    // print the premium value
    std::cout << "********************************************************"<<std::endl;
    std::cout << "THE SIMULATION DETAILS : "<<std::endl;
    std::cout << "The payoffs mean: "<< moypayoff << std::endl;
    std::cout << "The premium of the call option is : "<<premium << std::endl;
    //std::cout << "The payoffs std_deviation : "<< pay_stddev << std::endl;
    std::cout << "confidence interval of the mean estimation: [" << premium-2*(pay_stddev/sqrt(M)) << " ; "<< premium+2*(pay_stddev/sqrt(M)) << "]" << std::endl;
    std::cout << "The confidence interval size: "<< (premium+2*(pay_stddev/sqrt(M))) - (premium-2*(pay_stddev/sqrt(M))) << std::endl;


    // ****************************************************************************************************************************************
    // *************************************************** REDUCTION DE VARIANCE **************************************************************

    // ***************************************************** 1 VARIABLE ANTITETHIC *************************************************************

    // 2 Payoffs simulations loop
    // In this loop we replace half of the payoff already simulated with the antithetic variable

    for(int j(0);j<M;j++){

        // 3 path simulation loop
        for(int i(0);i<N+1;i++){
        S[i+1] = S[i]*exp((r-(0.5*sigma*sigma))*dt+sigma*sqrt(dt)*-alea_normal_dist(stddev)); // sim with antitethic variate -X of X
        }

        // 4 compute the sum of payoffs of the simulations
        payoffs[j] =  (payoffs[j] + std::max(S[N] - K,0.0))/2;
    }

    //double moypayoff_ant(moyenne<double,double>(payoffs,M)); // simulated payoff mean estimation with antithetic variates
    double moypayoff_rant((moyenne<double,double>(payoffs,M))); // moyenne des payoffs obtenus par variable antithétiques
    premium = exp(-r*T)*(moypayoff_rant); // Actualise la moyenne des payoff pour fournir le prix
    double pay_stddev_ant(0); // standard deviation of the payoffs simulated with antithetic variates
    pay_stddev_ant = sqrt(variance<double,double>(payoffs,M)); // standard deviation of the payoffs simulated

    // print the premium value
    std::cout << "********************************************************"<<std::endl;
    std::cout << "THE ANTITHETIC VARIATE OPTIMISATION SIMULATION DETAILS : "<<std::endl;
    std::cout << "The payoffs mean: "<< moypayoff_rant << std::endl;
    std::cout << "The premium of the call option is : "<<premium << std::endl;
    //std::cout << "The payoffs std_deviation : "<< pay_stddev_ant << std::endl;
    std::cout << "confidence interval of the mean estimation: [" << premium-2*(pay_stddev_ant/sqrt(M)) << " ; "<< premium+2*(pay_stddev_ant/sqrt(M)) << "]" << std::endl;
    std::cout << "The confidence interval size: "<< (premium+2*(pay_stddev_ant/sqrt(M)))-(premium-2*(pay_stddev_ant/sqrt(M))) << std::endl;


    // ****************************************************************************************************************************************
    // ******************************************************** 2 CONTROL VARIATE TECHNIQUE ******************************************************
    // Simulation of put payoffs to valuate the call

    for(int j(0);j<M;j++){

        // 3 path simulation loop
        for(int i(0);i<N+1;i++){
        S[i+1] = S[i]*exp((r-(0.5*sigma*sigma))*dt+sigma*sqrt(dt)*alea_normal_dist(stddev));
        }

        // 4 compute the sum of payoffs of the simulations
        payoffs[j] = exp(r*T)*S0 - K + std::max(K - S[N] ,0.0); // PUT PAYOFF
    }

    double callpayoff(0);
    callpayoff = moyenne<double,double>(payoffs,M); // CALL VALUE DEDUCES FROM CALL-PUTT PARITY FORMULA
    premium = exp(-r*T)*(callpayoff); // Actualise la moyenne des payoff pour fournir le prix
    double pay_stddev_varcontr(0); // standard deviation of the payoffs simulated with antithetic variates
    pay_stddev_varcontr = sqrt(variance<double,double>(payoffs,M)); // standard deviation of the payoffs simulated

    std::cout << "********************************************************"<<std::endl;
    std::cout << "THE CONTROL VARIATE OPTIMISATION DETAILS : "<<std::endl;
    std::cout << "The payoffs mean: "<< callpayoff << std::endl;
    std::cout << "The premium of the call option is : "<< premium << std::endl;
    //std::cout << "The payoffs std_deviation : "<< pay_stddev_varcontr << std::endl;
    std::cout << "confidence interval of the mean estimation: [" << premium-2*(pay_stddev_varcontr/sqrt(M)) << " ; "<< premium+2*(pay_stddev_varcontr/sqrt(M)) << "]" << std::endl;
    std::cout << "The confidence interval size: "<< (premium+2*(pay_stddev_varcontr/sqrt(M)))-(premium-2*(pay_stddev_varcontr/sqrt(M))) << std::endl;

    //verifying that the program finish
    return 0; //verifying that the program finish
}

// fonctions code

// Generate a number with follow a uniform distribution
// : http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
double alea_uniform_dist(){
    // generate random uniform numbers by c++11
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    //Use dis to transform the random unsigned int generated by gen into a double in [0, 1)
    return dis(gen); //Each call to dis(gen) generates a new random double
}


// Generate a number with follow a gaussian distribution
// : http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
double alea_normal_dist(const double& stddev){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dis(0.0, stddev);
    return dis(gen);
}
