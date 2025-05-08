#ifndef SPATIALSSA_MODEL_HPP
#define SPATIALSSA_MODEL_HPP

#define _USE_MATH_DEFINES

#define MAX_NEIGHBOURS    7
#define MAX_SPECIES       10
#define MAX_CELLS         10000

#define MAX_REACTIONS     20
#define MAX_REACTION_SIZE 5

const int MAX_EVENTS = MAX_CELLS*(MAX_REACTIONS + MAX_NEIGHBOURS*MAX_SPECIES);
// const int MAX_EVENTS = MAX_DIM*MAX_CELLS_PER_DIM*(MAX_REACTIONS + MAX_NEIGHBOURS*MAX_SPECIES);

#include <string>
#include <vector>
#include <map>          // For easily finding distance between cells if you only have absolute neighbour index
#include <set>          // For updating only some cells
#include <array>

#include <iostream>     // For writing to and reading from console, cout, cin, endl
#include <iomanip>      // For setting field width and alignment 

#include <ios>          // For specifying write type in streams
#include <fstream>      // For opening files
#include <sstream>      // For parsing strings from input files

#include <algorithm>

#include <chrono>       // For timing parts of code

#include <cmath>
#include <limits>


using std::string;
using std::to_string;
using std::stod;
using std::getline;

using std::vector;
using std::map;
using std::set;
using std::array;

// iostream
using std::cout;
using std::cin;
using std::endl;

// iomanip
using std::setw;
using std::left;

// ios
using std::ios;

// fstream
using std::ofstream;
using std::ifstream;

// sstream
using std::stringstream;

// algorithm
using std::min;
using std::max;

// limits
using std::numeric_limits;


// We make a Reaction class with the species index, state change, number of reactants, reaction name etc.
// Hence we have a Reaction for each actual reaction, and one corresponding to the diffusion of each species.
// We have an Event class containing only target cells and an index or pointer referencing the original Reaction that it refers to.


// Generic reaction
class Reaction {

    public:

        // Fixed properties of the reaction
        string reaction_name;
        bool is_diffusion_like;
        int reaction_species[MAX_REACTION_SIZE];
        int state_change[MAX_REACTION_SIZE];

        // Total number of nonunique species either created or destroyed. I.e. number of state changes that must be executed
        int n_reaction_species;

        // The species indices that react and determine the propensity. Index is the index in species list in reactants.
        int reactants[MAX_REACTION_SIZE];
        int n_reactants;

        // Rate index. Refers to rates stored in model
        int k_idx;

        // Number of events of this reaction that have occured
        int count;

        // Constructor
        Reaction();
};



class Event {

    public:

        // The position in the propensities list
        int idx;

        int target_cells[MAX_REACTION_SIZE];
        double base_rate;
        double rate;
        
        Reaction* reaction_ptr;

        // 0 is diffusion, 1 is reaction.
        
};


struct Vertex {

    double x, y, z;
    Vertex(double x, double y, double z);
    Vertex();

    Vertex operator+(const Vertex &v2);
    Vertex operator-(const Vertex &v2);
    Vertex operator/(const float &div);
    Vertex& operator/=(const float &div);

    Vertex cross(const Vertex &v2);
    double dot(const Vertex &v2);

    double distance(const Vertex &neighbour);
    double abs();
};

struct Face;


class Cell {

    public:

        // ID
        int idx;

        // Position of centre of the cell
        Vertex pos;

        // Position of edges
        vector<Vertex> vertices;
        
        // Idxs of ordered neighbours, forming a circle.
        vector<int> ordered_neighbours;

        // Faces in which the cell is a vertex.
        vector<Face*> adjacent_faces;

        // Neighbours
        vector<int> neighbours;
        map<int, double> neighbour_dist;
        map<int, double> neighbour_interface;
        unsigned int n_neighbours;
        
        // Current number of each species in the cell
        int state[MAX_SPECIES];

        // The events that are affected by a change of species in this cell
        vector<Event*> affected_events;
        unsigned int n_affected_events;

        // Volume and area for correcting the propensities
        double size;

        // The number of particles that it can hold and the number of particles it holds now.
        double capacity;
        double content;

        double phi;

        // The dimension of the cell (is it a surface or a volume)
        int dim;

        Cell();
        void set_species_number(int species_idx, int value);
        void add_neighbour(int neighbour_idx, double distance, double interface);
        void add_neighbour_alt_dim(int neighbour_idx, double interface);

        bool has_neighbour(int neighbour_idx);

        void add_affected_event(Event* event_ptr);

        double distance(const Cell &neighbour);

};

struct Face {

    vector<Cell*> cells;
    // vector<Vertex> edge_normals;

    // Only works if faces are triangles
    // For each vertex, we have a theta term l1*l2*cos(theta_3).
    // These are positive for triangles, since their edges form sharp angles.
    vector<double> theta_terms;
    vector<double> lengths;
    double area;

    Face(Cell *c1, Cell *c2, Cell *c3);
    // void generate_edge_normals();
    void generate_angles_and_lengths();

};

struct Edge{
    Cell* c1;
    Cell* c2;
    Edge(Cell &c1, Cell &c2);
};

class Model {

    public:

        // Safety check parameters:
        double F_max;
        int n_reject;
        int n_recalculate;

        // LLPS parameters for protein
        // double particles_per_area = 80000;
        // double particles_per_volume = 30000000;

        // For phase seperation
        double particles_per_area;
        double particles_per_volume;

        // For crowding. The fraction of volume or membrane area in the cell that is unoccupied by crowders
        // Affects only normal reactions and boundary reactions. Not diffusion.
        double free_volume_fraction;
        double free_area_fraction;
        map<int, double> free_size_fraction;
        

        // The Flory-Huggins chi that determines whether an interaction is favourable or unfavourable.
        double chi;

        // The strength of the interfacial energy caused by surface tension.
        double kappa;


        // Name of the project i.e. run, etc.
        string dir_name;
        
        string run_name;
        string project_name;
        string project_date;

        // Name of the model i.e. Lotka-Volterra
        string model_name;

        // Number of repeats
        int n_repeats;

        // Number of the run (for parameter sweeps)
        int i_run;

        // Apply phase seperation effects
        bool phase_separate;

        // Apply crowding effects
        bool crowding;


        // Bool to check if initial recording has already happened
        bool initial_recording;

        bool record_indv_cells;
        bool record_indv_species;
        bool record_reactions;
        bool record_propensities;

        // Bool for verbose printing of propensities
        bool verbose;
        bool timing;

        // Time and iteration
        double t;
        int i_iter;

        double t_max;
        long long int n_iterations, break_iter;

        // The cells that make up the grid
        // map<int, array<Cell, MAX_CELLS_PER_DIM>> cell_map;
        Cell cells[MAX_CELLS];
        int n_cells;
        map<int, int> n_cells_per_dim;

        // All the faces of the mesh
        vector<Face> faces;

        // Normal reactions that can occur, one per reaction mechanism.
        Reaction reactions[MAX_REACTIONS]; 
        int n_reactions;

        // Reactions that cross dimensionality boundaries, such as membrane attachment.
        Reaction boundary_reactions[MAX_REACTIONS]; 
        int n_boundary_reactions;

        // Diffusion reactions that can occur, one per species.
        Reaction diffusion_reactions[MAX_SPECIES];
        // Number is counted by n_species

        
        // The rate constants of the possible reactions
        map<string, int> rate_constant_map;
        string rate_constant_names[MAX_REACTIONS + MAX_SPECIES];
        double rate_constants[MAX_REACTIONS + MAX_SPECIES];
        int n_rate_constants;

        int initial_species_number[MAX_SPECIES];

        // The chemical species, reactants and products.
        string species[MAX_SPECIES];
        int species_dimensions[MAX_SPECIES];
        int n_species;

        // All possible events and their propensities, i.e. all possible combinations of cells, reactions, diffusions and directions.
        Event events[MAX_EVENTS];
        int event_count[MAX_REACTIONS + MAX_SPECIES];
        int n_events;
        double propensities[MAX_EVENTS];
        double a0;

        // The last performed reaction and random number
        int mu;
        vector<int> mu_history;
        double r_mu;

        // Vectors containing values to sweep over
        // The elements of each vector are vectors which contain a range of species numbers for one species
        // The index vector contains the species they refer to. So the first dimensons are equal.

        // Vectors for sweep over the abundance of species
        vector<vector<int>> sweeps_species_number;
        vector<int> sweeps_species_idx;

        // Vectors for sweep over the rate coefficients and diffusion constants
        vector<vector<double>> sweeps_rate_constants;
        vector<int> sweeps_rate_constant_idx;

        
        // Constructor
        Model();

        // TODO: create new function which assigns reactants. Now there is code duplication in add_species and add_reaction.

        void create_cartesian_grid(const double d_x, const double d_y, const double d_z, const int n_x, const int n_y, const bool periodic);
        void create_spherical_grid(const int n_points, const double r);

        void run();
        void do_parameter_sweep();

        void load_model(const string file_name);
        void load_initial_conditions(const string file);

        void define_diffusion_reactions();
        void define_normal_reactions();
        void define_boundary_reactions();

        void record_concentrations(const int n_iter);
        void record_mu(const int n_iter);
        void record_parameters();
        void record_initial_conditions();
        void record_event_data();
        void record_cell_data();
        void error_dump();

        void record_all_species_in_cell(const int n_iter, const int i_cell);
        void record_all_cells_for_species(const int i_species, const int n_iter);
        void record_all_reactions(const int n_iter);
        void record_all_propensities(const int n_iter);


        // Only public because it is needed by load_initial conditions
        int get_species_idx(const string species_name);
        int get_rate_idx(const string rate_constant_name);
        void distribute_particles(const int n_particles, int const i_species);
        void make_particle_spot(const int cell_idx, const double spot_radius, const int count, const int i_species);

        void add_species_sweep(const int i_species, const int n_start, const int n_end, const int n_points, const string sweep_type);
        void add_rate_constant_sweep(const int i_rate_constant, const double k_start, const double k_end, const int n_points, const string sweep_type);

        void load_phase_separation_parameters();


    private:

        void add_species(const string species_name, const int dimension, const string dif_coefficient_name);
        void add_reaction(const string reaction_species[], const int state_changes[], const string rate_constant_name, const int n_reaction_species, const string reaction_name);
        void add_rate_constant(const string rate_constant_name);

        void reset();

        void calculate_propensities();
        void recalculate_propensities();

        double get_event_propensity(const Event* event_ptr);

        //LLPS
        void update_diffusion_rate(Event* event_ptr, const double del_f);


        // TEMP
        // Calculating total gradient squared energy
        // double calculate_gradient_energy();
        // Calculating the difference of a gradient squared energy. Interface like calculate_energy_dif.
        // double calculate_gradient_energy_dif(const Event* event_ptr);
        // End TEMP

        double calculate_energy_dif(const Event* event_ptr);
        double calculate_enthalpy_dif(const double phi, const double del_phi);
        double calculate_enthalpy_dif_adapted(const double phi_o, const double phi_t);
        // double calculate_surface_tension_dif_grid(const double phi, const double del_phi, Cell *cell_ptr);
        double calculate_gradient_energy_dif_mesh(const double phi_o, const double del_phi_o, const double del_phi_t, Cell *cell_ptr_o, Cell *cell_ptr_t);
        double calculate_repulsion_dif(const double phi, const double del_phi);
        // End LLPS

        void update_state();

};





#endif //SPATIALSSA_MODEL_HPP