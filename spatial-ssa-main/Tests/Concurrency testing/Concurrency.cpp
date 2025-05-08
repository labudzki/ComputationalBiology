#include <iostream>
#include <vector>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;

using std::vector;

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;


#include <random>

using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;

// Random number generation
// Generation
random_device rd;  							// Will be used to obtain a seed for the random number engine. As rd() not guaranteed to follow any distribution.
mt19937 gen(rd());							// Standard mersenne_twister_engine seeded with rd(). Is statistically random.
uniform_real_distribution<> dis(0.0, 1.0); 	// Call to "dis(gen)" will now generate a random double in [0,1)


int main(){

    high_resolution_clock sum_clock;
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();

    // Getting number of milliseconds as an integer.
    auto calculate_time = milliseconds(0);



    int n_props = 50*50 * (4*10 + 10);
    int n_repeats = 100000;
    // int n_props = 10;
    double r, sum;
    auto propensities = new double[n_props];

    // Generating array
    for(int i_prop = 0; i_prop < n_props; i_prop++){
        r = dis(gen);
        propensities[i_prop] = r;

        // cout << propensities[i_prop] << endl;
        // cin.get();
    }

    // Summing over array
    for (int i_repeat = 0; i_repeat < n_repeats; i_repeat++){
        if (i_repeat%1000 == 0){cout << i_repeat/1000 << endl;}
        sum = 0;

        // Starting clock
        t1 = sum_clock.now();
        for(int i_prop = 0; i_prop < n_props; i_prop++){
            sum += propensities[i_prop];
        }

        // Ending clock
        t2 = sum_clock.now();
        calculate_time += duration_cast<milliseconds>(t2 - t1);
        
    }

    cout << "Time total: " << calculate_time.count() << endl;


    delete[] propensities;

    cout << "hi" << endl;
    cin.get();

    return 0;
}