#include <iostream>
using std::cout, std::endl, std::cin;

#include <random>
using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;

// Random number generation
// Generation
random_device rd;  							// Will be used to obtain a seed for the random number engine. As rd() not guaranteed to follow any distribution.
mt19937 gen(rd());							// Standard mersenne_twister_engine seeded with rd(). Is statistically random.
uniform_real_distribution<> dis(0.0, 1.0); 	// Call to "dis(gen)" will now generate a random double in [0,1)

#include <string>
using std::string, std::to_string;

#include <fstream>
using std::ofstream;

#include <ios>
using std::ios;

#include <string>
using std::string;


 
int main(int argc, char** argv)
{   

    // Loading in directory name
    string dir_name = "";
    if(argc > 1){
        dir_name = argv[1];
    } else{
        dir_name = "Data";
    }

    string project_date, project_name;
    // Setting project date. Time the model is run
    time_t start_time = time(0);
    tm *ltm = localtime(&start_time);

    project_date = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "_"
                 + to_string(ltm->tm_hour) + "-" + to_string(ltm->tm_min) + "_" + to_string(ltm->tm_sec);

    project_name = project_date;

    cout << "Project name: " << project_name << endl;

    string file_name, file_path;
    file_name = project_name + "_numbers.txt";
    file_path = dir_name + "/" + file_name;

    cout << "Directory name: " << dir_name << endl;
    cout << "File path: " << file_path << endl;
    cout << "Random device: " << rd() << endl;
    cout << "Random number: " << dis(gen) << endl;

	ofstream datafile;
	
	datafile.open(file_path, ios::app);
    for(int i = 0; i<10; i++){
        datafile << dis(gen) << "\n";
    }

	datafile << "\n";
	datafile.close();

    // fs::current_path(fs::temp_directory_path());
    // fs::create_directories("sandbox/1/2/a");
    // fs::create_directory("sandbox/1/2/b");
    // fs::permissions("sandbox/1/2/b", fs::perms::others_all, fs::perm_options::remove);
    // fs::create_directory("sandbox/1/2/c", "sandbox/1/2/b");
    // // std::system("ls -l sandbox/1/2");
    // // std::system("tree sandbox");
    // fs::remove_all("sandbox");

}